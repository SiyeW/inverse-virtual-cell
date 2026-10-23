/* Scoped compatibility for an already-satisfied fchmod on CIFS.
 * Does not grant permissions or ignore a requested permission change.
 * Enable only in an independently owned test/run with an explicit output root.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/vfs.h>
#include <unistd.h>

int fchmod(int fd, mode_t requested) {
    int (*real_fchmod)(int, mode_t) = dlsym(RTLD_NEXT, "fchmod");
    if (!real_fchmod) { errno = ENOSYS; return -1; }
    int result = real_fchmod(fd, requested);
    int saved_errno = errno;
    if (result == 0 || (saved_errno != EPERM && saved_errno != EACCES)) return result;

    const char *configured = getenv("MRO_CIFS_SAME_MODE_ROOT");
    char root[PATH_MAX], proc[64], target[PATH_MAX];
    struct stat st;
    struct statfs fs;
    if (!configured || !realpath(configured, root) || strlen(root) < 2) goto denied;
    if (fstat(fd, &st) != 0 || fstatfs(fd, &fs) != 0) goto denied;
    if (!S_ISREG(st.st_mode) && !S_ISDIR(st.st_mode)) goto denied;
    if ((st.st_mode & 07777) != (requested & 07777)) goto denied;
    if ((unsigned long)fs.f_type != 0xff534d42UL &&
        (unsigned long)fs.f_type != 0xfe534d42UL) goto denied;
    snprintf(proc, sizeof(proc), "/proc/self/fd/%d", fd);
    ssize_t n = readlink(proc, target, sizeof(target) - 1);
    if (n < 0 || n >= (ssize_t)sizeof(target) - 1) goto denied;
    target[n] = '\0';
    size_t length = strlen(root);
    if (strncmp(target, root, length) || target[length] != '/') goto denied;
    /* Descriptor identity, filesystem, directory scope and requested state have
       all been checked. No metadata or permission is changed by this branch. */
    fprintf(stderr, "[cifs-same-mode] already %04o: %s\n", requested & 07777, target);
    errno = saved_errno;
    return 0;
denied:
    errno = saved_errno;
    return result;
}
