import htsjdk.samtools.*;
import com.google.gson.*;
import java.nio.file.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.io.*;
import java.util.zip.GZIPOutputStream;
import java.security.MessageDigest;

/** Read-only indexed BAM access through HTSJDK; never creates a BAM. */
public class GuidedRnaReadCache {
    static final Gson JSON = new GsonBuilder().disableHtmlEscaping().create();
    static String sha(Path path) throws Exception {
        MessageDigest md=MessageDigest.getInstance("SHA-256");
        try(InputStream in=Files.newInputStream(path)) {
            byte[] b=new byte[65536];int n;
            while((n=in.read(b))!=-1)md.update(b,0,n);
        }
        return HexFormat.of().formatHex(md.digest());
    }
    static void atomic(Path path, JsonObject object) throws Exception {
        Path temp=path.resolveSibling(path.getFileName()+".tmp");
        Files.write(temp,JSON.toJson(object).getBytes(StandardCharsets.UTF_8));
        Files.move(temp,path,StandardCopyOption.REPLACE_EXISTING,StandardCopyOption.ATOMIC_MOVE);
    }
    static void verifyFile(String name, JsonObject identity) throws Exception {
        Path p=Paths.get(name);
        if(Files.size(p)!=identity.get("bytes").getAsLong() ||
           Files.getLastModifiedTime(p).to(TimeUnit.NANOSECONDS)!=identity.get("mtime_ns").getAsLong())
            throw new IllegalStateException("Source identity changed: "+name);
    }
    static String sample(SAMRecord read, JsonObject region, Map<String,String> calls, Set<String> ntc) {
        String cb=read.getStringAttribute("CB");if(cb==null)return null;
        String guide=calls.get(cb);if(guide==null)return null;
        if(ntc.contains(guide))return "NTC";
        for(JsonElement element:region.getAsJsonArray("guides")) {
            JsonObject g=element.getAsJsonObject();
            if(guide.equals(g.get("guide").getAsString()) &&
               read.getAlignmentStart()-1<g.get("end0").getAsInt()+100 &&
               read.getAlignmentEnd()>Math.max(0,g.get("start0").getAsInt()-100))
                return g.get("sample_id").getAsString();
        }
        return null;
    }
    public static void main(String[] args) throws Exception {
        JsonObject cfg=JsonParser.parseString(Files.readString(Paths.get(args[0]))).getAsJsonObject();
        verifyFile(cfg.get("bam").getAsString(),cfg.getAsJsonObject("bam_identity"));
        verifyFile(cfg.get("bam_index").getAsString(),cfg.getAsJsonObject("bam_index_identity"));
        SamInputResource input=SamInputResource.of(Paths.get(cfg.get("bam").getAsString()).toFile())
            .index(Paths.get(cfg.get("bam_index").getAsString()).toFile());
        try(SamReader reader=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(input)) {
            JsonArray expected=cfg.getAsJsonArray("reference_contigs_ordered");
            List<SAMSequenceRecord> actual=reader.getFileHeader().getSequenceDictionary().getSequences();
            if(actual.size()!=expected.size())throw new IllegalStateException("Reference cardinality differs");
            for(int i=0;i<actual.size();i++) {
                JsonObject e=expected.get(i).getAsJsonObject();
                if(!actual.get(i).getSequenceName().equals(e.get("name").getAsString()) ||
                   actual.get(i).getSequenceLength()!=e.get("length").getAsInt())
                    throw new IllegalStateException("Complete ordered reference identity differs");
            }
            if(cfg.get("mode").getAsString().equals("probe")) {
                JsonObject region=cfg.getAsJsonArray("regions").get(0).getAsJsonObject();
                int n=0,cb=0,ub=0;
                try(SAMRecordIterator it=reader.queryOverlapping(region.get("chrom").getAsString(),region.get("start0").getAsInt()+1,region.get("end0").getAsInt())) {
                    while(it.hasNext() && n<20){SAMRecord r=it.next();n++;if(r.getAttribute("CB")!=null)cb++;if(r.getAttribute("UB")!=null)ub++;}
                }
                JsonObject result=new JsonObject();result.addProperty("records",n);result.addProperty("CB_present",cb);result.addProperty("UB_present",ub);result.addProperty("reference_inventory_matches",true);System.out.println(JSON.toJson(result));return;
            }
            Map<String,String> calls=new HashMap<>();
            for(String line:Files.readAllLines(Paths.get(cfg.get("assignments_tsv").getAsString()),StandardCharsets.UTF_8)) {
                String[] a=line.split("\t",-1);if(a.length!=2||calls.put(a[0],a[1])!=null)throw new IllegalStateException("Ambiguous barcode assignment");
            }
            Set<String> ntc=new HashSet<>(Files.readAllLines(Paths.get(cfg.get("NTC_guides_txt").getAsString()),StandardCharsets.UTF_8));
            Path root=Paths.get(cfg.get("local_staging").getAsString());Files.createDirectories(root);
            JsonObject identity=cfg.getAsJsonObject("identity");
            Set<String> completed=new HashSet<>();
            for(JsonElement x:cfg.getAsJsonArray("skip_regions"))completed.add(x.getAsString());
            int done=completed.size();
            long cacheBytes=0;
            for(String rid:completed)cacheBytes+=Files.size(root.resolve(rid+".jsonl.gz"));
            boolean waitForUpload=cfg.has("wait_for_upload") && cfg.get("wait_for_upload").getAsBoolean();
            for(JsonElement element:cfg.getAsJsonArray("regions")) {
                JsonObject region=element.getAsJsonObject();String rid=region.get("region_id").getAsString();
                if(completed.contains(rid))continue;
                Path file=root.resolve(rid+".jsonl.gz"),partial=root.resolve(rid+".partial.gz");
                if(Files.exists(file)||Files.exists(partial))throw new IllegalStateException("Unaccepted artifact exists: "+rid);
                long start=System.nanoTime();int fetched=0,selected=0,missing=0;boolean capped=false;
                JsonObject status=new JsonObject();status.addProperty("stage","extract_region");status.addProperty("region",rid);status.addProperty("completed_regions",done);atomic(root.resolve("source_status.json"),status);
                try(BufferedWriter out=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(Files.newOutputStream(partial)),StandardCharsets.UTF_8));
                    SAMRecordIterator it=reader.queryOverlapping(region.get("chrom").getAsString(),region.get("start0").getAsInt()+1,region.get("end0").getAsInt())) {
                    JsonObject header=new JsonObject();header.addProperty("kind","RNA_SAM_cache_header");header.add("identity",identity);header.add("region",region);out.write(JSON.toJson(header));out.newLine();
                    while(it.hasNext()) {
                        SAMRecord r=it.next();fetched++;
                        if(fetched>cfg.get("max_fetched_records_per_region").getAsInt()){capped=true;break;}
                        if(r.getReadUnmappedFlag()||r.isSecondaryAlignment()||r.getSupplementaryAlignmentFlag()||r.getReadFailsVendorQualityCheckFlag()||r.getMappingQuality()<30)continue;
                        String label=sample(r,region,calls,ntc);if(label==null)continue;
                        selected++;if(r.getAttribute("UB")==null)missing++;
                        JsonObject row=new JsonObject();row.addProperty("kind","original_SAM_record");row.addProperty("sample",label);row.addProperty("sam",r.getSAMString().stripTrailing());out.write(JSON.toJson(row));out.newLine();
                    }
                }
                Files.move(partial,file,StandardCopyOption.ATOMIC_MOVE);
                JsonObject ready=new JsonObject();ready.addProperty("stage","region_ready_for_verified_upload");ready.add("region",region);ready.addProperty("bytes",Files.size(file));ready.addProperty("sha256",sha(file));ready.addProperty("fetched_records",fetched);ready.addProperty("selected_records",selected);ready.addProperty("selected_without_UB",missing);ready.addProperty("scan_capped",capped);ready.addProperty("seconds",(System.nanoTime()-start)/1e9);atomic(root.resolve(rid+".ready.json"),ready);
                cacheBytes+=Files.size(file);
                if(cacheBytes>cfg.get("max_local_staging_bytes").getAsLong())throw new IllegalStateException("Local cache safety bound exceeded; ready regions remain preserved");
                // Backpressure keeps only one small region cache on this host.
                while(waitForUpload && !Files.exists(root.resolve(rid+".uploaded"))) {
                    if(Files.exists(root.resolve("abort")))throw new IllegalStateException("Upload coordinator stopped; preserve the ready cache");
                    Thread.sleep(200);
                }
                done++;status.addProperty("stage",waitForUpload?"region_uploaded":"region_cached");status.addProperty("completed_regions",done);atomic(root.resolve("source_status.json"),status);
            }
            verifyFile(cfg.get("bam").getAsString(),cfg.getAsJsonObject("bam_identity"));
            verifyFile(cfg.get("bam_index").getAsString(),cfg.getAsJsonObject("bam_index_identity"));
            JsonObject finalStatus=new JsonObject();finalStatus.addProperty("stage","source_cache_complete");finalStatus.addProperty("completed_regions",done);atomic(root.resolve("source_status.json"),finalStatus);
        }
    }
}
