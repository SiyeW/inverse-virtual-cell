//go:build linux

// Opt-in CIFS log identity proof. No scientific input or output is changed.
package util

import (
	"bytes"
	"crypto/rand"
	"encoding/hex"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"syscall"
)

const cifsLogMagic = 0xff534d42
const smb2LogMagic = 0xfe534d42

// Only a positively identified CIFS descriptor may use the alternate proof.
// Native files and non-opted-in runs retain upstream's SameFile requirement.
func verifyCifsLogIdentity(f *os.File, previous, current os.FileInfo) error {
	if os.Getenv("MRO_CIFS_COMPAT") != "1" {
		return fmt.Errorf("CIFS log identity proof is not enabled")
	}
	var fs syscall.Statfs_t
	if err := syscall.Fstatfs(int(f.Fd()), &fs); err != nil {
		return err
	}
	if uint64(fs.Type) != cifsLogMagic && uint64(fs.Type) != smb2LogMagic {
		return fmt.Errorf("log descriptor is not on CIFS/SMB2 (type %#x)", fs.Type)
	}
	return proveCifsLogPath(f, previous, current)
}

// The caller holds the logger mutex. Appending an unpredictable marker through
// the original descriptor, then reading it at the same offset through the path,
// proves live shared content. A copied/replaced log cannot acquire this marker.
// A metadata-only, size-only or matching-prefix check would not prove this.
func proveCifsLogPath(f *os.File, previous, current os.FileInfo) error {
	root := os.Getenv("MRO_CIFS_SAME_MODE_ROOT")
	if root == "" || !filepath.IsAbs(root) || !filepath.IsAbs(f.Name()) || filepath.Base(f.Name()) != "_log" {
		return fmt.Errorf("explicit absolute owned log root required")
	}
	root, err := filepath.EvalSymlinks(root)
	if err != nil {
		return err
	}
	path, err := filepath.EvalSymlinks(f.Name())
	if err != nil {
		return err
	}
	rel, err := filepath.Rel(root, path)
	if err != nil || rel == ".." || strings.HasPrefix(rel, ".."+string(os.PathSeparator)) {
		return fmt.Errorf("log is outside the explicit owned root")
	}
	linkInfo, err := os.Lstat(f.Name())
	if err != nil {
		return err
	}
	if !previous.Mode().IsRegular() || !current.Mode().IsRegular() || !linkInfo.Mode().IsRegular() {
		return fmt.Errorf("regular non-symlink log required")
	}
	oldStat, okOld := previous.Sys().(*syscall.Stat_t)
	newStat, okNew := current.Sys().(*syscall.Stat_t)
	if !okOld || !okNew || oldStat.Dev != newStat.Dev {
		return fmt.Errorf("log filesystem changed")
	}
	var nonce [32]byte
	if _, err := rand.Read(nonce[:]); err != nil {
		return err
	}
	marker := []byte("\n[cifs-log-identity] " + hex.EncodeToString(nonce[:]) + "\n")
	if n, err := f.Write(marker); err != nil {
		return err
	} else if n != len(marker) {
		return io.ErrShortWrite
	}
	// The descriptor offset after O_APPEND is authoritative even when CIFS
	// returns a stale size through stat. The logger mutex excludes local writes.
	end, err := f.Seek(0, io.SeekCurrent)
	if err != nil || end < int64(len(marker)) {
		return fmt.Errorf("cannot locate appended log identity marker: %v", err)
	}
	offset := end - int64(len(marker))
	if err := f.Sync(); err != nil {
		return err
	}
	reader, err := os.Open(f.Name())
	if err != nil {
		return err
	}
	defer reader.Close()
	actual := make([]byte, len(marker))
	if _, err := reader.ReadAt(actual, offset); err != nil {
		return fmt.Errorf("original descriptor marker not visible through log path: %w", err)
	}
	if !bytes.Equal(actual, marker) {
		return fmt.Errorf("original descriptor and log path have different live content")
	}
	return nil
}
