//go:build !linux

package util

import (
	"fmt"
	"os"
)

func verifyCifsLogIdentity(f *os.File, previous, current os.FileInfo) error {
	return fmt.Errorf("CIFS log identity proof is supported only on Linux")
}
