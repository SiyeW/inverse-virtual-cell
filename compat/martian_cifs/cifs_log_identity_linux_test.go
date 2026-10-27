//go:build linux

package util

import (
	"os"
	"path/filepath"
	"strings"
	"sync"
	"syscall"
	"testing"
)

type changedCifsInode struct {
	os.FileInfo
	stat syscall.Stat_t
}

func (i changedCifsInode) Sys() interface{} { return &i.stat }

func changedLogInfo(t *testing.T, f *os.File) os.FileInfo {
	t.Helper()
	i, err := f.Stat()
	if err != nil {
		t.Fatal(err)
	}
	s := *i.Sys().(*syscall.Stat_t)
	s.Ino++
	return changedCifsInode{FileInfo: i, stat: s}
}

func logProofFixture(t *testing.T, root string) (*os.File, string) {
	t.Helper()
	t.Setenv("MRO_CIFS_SAME_MODE_ROOT", root)
	p := filepath.Join(root, "_log")
	f, err := os.OpenFile(p, os.O_WRONLY|os.O_CREATE|os.O_APPEND, 0600)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { f.Close() })
	if _, err := f.WriteString("same initial log contents\n"); err != nil {
		t.Fatal(err)
	}
	return f, p
}

func TestCifsLogProofRemappedIdentity(t *testing.T) {
	f, p := logProofFixture(t, t.TempDir())
	previous := changedLogInfo(t, f)
	current, _ := os.Stat(p)
	if os.SameFile(previous, current) {
		t.Fatal("test must simulate unequal file IDs")
	}
	if err := proveCifsLogPath(f, previous, current); err != nil {
		t.Fatal(err)
	}
	b, _ := os.ReadFile(p)
	if !strings.Contains(string(b), "[cifs-log-identity]") {
		t.Fatal("live-content proof was not written")
	}
}

func TestCifsLogProofRejectsIdenticalReplacement(t *testing.T) {
	root := t.TempDir()
	f, p := logProofFixture(t, root)
	previous, _ := f.Stat()
	old, _ := os.ReadFile(p)
	if err := os.Rename(p, filepath.Join(root, "old_log")); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(p, old, 0600); err != nil {
		t.Fatal(err)
	}
	current, _ := os.Stat(p)
	if err := proveCifsLogPath(f, previous, current); err == nil {
		t.Fatal("copied log replacement was accepted")
	}
}

func TestCifsLogProofRejectsSymlink(t *testing.T) {
	root := t.TempDir()
	f, p := logProofFixture(t, root)
	previous, _ := f.Stat()
	moved := filepath.Join(root, "actual")
	if err := os.Rename(p, moved); err != nil {
		t.Fatal(err)
	}
	if err := os.Symlink(moved, p); err != nil {
		t.Fatal(err)
	}
	current, _ := os.Stat(p)
	if err := proveCifsLogPath(f, previous, current); err == nil {
		t.Fatal("symlink was accepted")
	}
}

func TestCifsLogProofRejectsOutsideRoot(t *testing.T) {
	f, p := logProofFixture(t, t.TempDir())
	t.Setenv("MRO_CIFS_SAME_MODE_ROOT", t.TempDir())
	i, _ := f.Stat()
	c, _ := os.Stat(p)
	if err := proveCifsLogPath(f, i, c); err == nil {
		t.Fatal("outside-root log was accepted")
	}
}

func TestVerifyLogRetainsStrictNativeChecks(t *testing.T) {
	root := t.TempDir()
	f, p := logProofFixture(t, root)
	before, _ := f.Stat()
	t.Setenv("MRO_CIFS_COMPAT", "1")
	current, _ := os.Stat(p)
	var fs syscall.Statfs_t
	if err := syscall.Fstatfs(int(f.Fd()), &fs); err != nil {
		t.Fatal(err)
	}
	if uint64(fs.Type) == cifsLogMagic || uint64(fs.Type) == smb2LogMagic {
		t.Skip("native test needs a non-CIFS tmpdir")
	}
	if err := verifyCifsLogIdentity(f, before, current); err == nil {
		t.Fatal("native filesystem used CIFS proof")
	}
	t.Setenv("MRO_CIFS_COMPAT", "")
	if err := verifyCifsLogIdentity(f, before, current); err == nil {
		t.Fatal("non-opted-in run used CIFS proof")
	}
}

func TestVerifyLogMissingAndReplaced(t *testing.T) {
	root := t.TempDir()
	f, p := logProofFixture(t, root)
	oldLogger := LOGGER
	LOGGER = &Logger{fileWriter: f}
	t.Cleanup(func() { LOGGER = oldLogger })
	if err := VerifyLogFile(); err != nil {
		t.Fatal(err)
	}
	if err := os.Remove(p); err != nil {
		t.Fatal(err)
	}
	if err := VerifyLogFile(); err == nil {
		t.Fatal("missing log was accepted")
	}
	if err := os.WriteFile(p, []byte("same initial log contents\n"), 0600); err != nil {
		t.Fatal(err)
	}
	if err := VerifyLogFile(); err == nil {
		t.Fatal("replacement log was accepted")
	}
}

func TestVerifyLogConcurrentWrites(t *testing.T) {
	f, _ := logProofFixture(t, t.TempDir())
	oldLogger := LOGGER
	LOGGER = &Logger{fileWriter: f}
	t.Cleanup(func() { LOGGER = oldLogger })
	var group sync.WaitGroup
	for i := 0; i < 4; i++ {
		group.Add(1)
		go func() {
			defer group.Done()
			for j := 0; j < 20; j++ {
				LOGGER.WriteString("concurrent entry\n")
			}
		}()
	}
	for j := 0; j < 20; j++ {
		if err := VerifyLogFile(); err != nil {
			t.Fatal(err)
		}
	}
	group.Wait()
}

func TestCifsLogActualStorage(t *testing.T) {
	parent := os.Getenv("MRO_CIFS_LOG_TEST_ROOT")
	if parent == "" {
		t.Skip("explicit disposable CIFS test root not supplied")
	}
	root, err := os.MkdirTemp(parent, "log-proof-")
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { os.RemoveAll(root) })
	f, p := logProofFixture(t, root)
	t.Setenv("MRO_CIFS_COMPAT", "1")
	current, _ := os.Stat(p)
	if err := verifyCifsLogIdentity(f, changedLogInfo(t, f), current); err != nil {
		t.Fatal(err)
	}
	previous, _ := f.Stat()
	old, _ := os.ReadFile(p)
	if err := os.Rename(p, filepath.Join(root, "preserved_old_log")); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(p, old, 0600); err != nil {
		t.Fatal(err)
	}
	current, _ = os.Stat(p)
	if err := verifyCifsLogIdentity(f, previous, current); err == nil {
		t.Fatal("actual CIFS replacement was accepted")
	}
}
