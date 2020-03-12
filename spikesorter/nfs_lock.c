#define _GNU_SOURCE
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <signal.h>
#include <string.h>
#include <errno.h>
#include <error.h>
#include <sys/time.h>
#include <time.h>
#include <sys/utsname.h>
#include <sys/file.h>


static char *
hostname (void)
{
  static struct utsname buf;
  uname (&buf);
  return buf.nodename;
}

static long unsigned
get_starttime (int pid)
{
  static char *line = NULL;
  static size_t len = 0;
  char *p;
  long unsigned starttime;
  char *statfilename;
  FILE *statfile;

  if (asprintf (&statfilename, "/proc/%d/stat", pid) == -1) exit (1);
  if ((statfile = fopen (statfilename, "r")) == NULL)
    return 0;
  if (getline (&line, &len, statfile) == -1) {
    fprintf (stderr, "Can't read %s\n", statfilename);
    exit (1);
  }
  p = strrchr (line, ')') + 2;
  if (sscanf (p, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %lu", &starttime) != 1)
    return 0;                   /* process must just have died */
  fclose (statfile);
  free (statfilename);
  return starttime;
}

long long
now (void)
{
  struct timeval t;
  gettimeofday (&t, 0); 
  return t.tv_sec * 1000000LL + t.tv_usec;
}

static char *filename;
static char *lockfile_name;

int
fail (char *fmt, ...)
{
  va_list ap;

  va_start (ap, fmt);
  vprintf (fmt, ap);
  va_end (ap);
  printf ("DISPATCH FAILED.\n");
  unlink (filename);
  exit (0);
}

int
main (int argc, char **argv)
{
  int fd;
  
  if (argc == 2 && strcmp (argv[1], "unique") == 0) {
    if (asprintf (&filename, "_%s_%d_%lu_%lld", hostname (), getppid (), get_starttime (getppid ()), now()) == -1) exit (1);
    printf ("%s\n", filename);
    return 0;
  }
  
  if (argc != 3)
    return 1;

  if (asprintf (&filename, "%s/%s", argv[1], argv[2]) == -1) exit (1);
  if (asprintf (&lockfile_name, "%s/lock_file", argv[1]) == -1) exit (1);
  fd = open (filename, O_RDWR | O_CREAT | O_TRUNC, 0644);
  if (fd == -1) {
    perror ("open");
    return 1;
  }
  if (write (fd, filename, strlen (filename)) == -1)
    error_at_line (1, errno, __FILE__, __LINE__, "write");
  close (fd);
  while (link (filename, lockfile_name) != 0) {
    struct stat s;
    stat (filename, &s);
    if (s.st_nlink == 1) {
      FILE *f;
      static char *line = NULL;
      static size_t len = 0;
      char *lf_hostname;
      pid_t lf_pid;
      long unsigned lf_starttime;
      long long lf_timeofday;
      int lfd;
      char *second_lock = "/var/local/lock/spikesort_control_panel";
      char *existing_lockfile;
      
      if ((lfd = open (second_lock, O_RDONLY | O_CREAT, 0644)) == -1)
        fail ("Can't open %s: %s\n", second_lock, strerror (errno));
      flock (lfd, LOCK_EX);

      if ((f = fopen (lockfile_name, "r")) == NULL)
        fail ("Can't open %s: %s\n", lockfile_name, strerror (errno));
      if (getline (&line, &len, f) == -1)
        fail ("Can't read %s\n", lockfile_name);
      fclose (f);
      existing_lockfile = strdup (line);
      {
        char *p;

        while ((p = strchr (line, '_')))
          *p = ' ';
        p = strrchr (line, '/');
        if (sscanf (p, "%*s %ms %d %ld %lld", &lf_hostname, &lf_pid, &lf_starttime, &lf_timeofday) != 4)
          fail ("Can't parse %s\n", lockfile_name);
      }
      if (strcmp (lf_hostname, hostname ()) == 0) {
        if (kill (lf_pid, 0) == 0 && lf_starttime == get_starttime (lf_pid))
          fail ("Process %d on host %s is dispatching this experiment.\n", lf_pid, lf_hostname);
        else {
          if (unlink (existing_lockfile) == -1)
            fail ("Can't unlink %s: %s\n", filename, strerror (errno));
          if (unlink (lockfile_name) == -1)
            fail ("Can't unlink %s: %s\n", lockfile_name, strerror (errno));
          flock (lfd, LOCK_UN);
          close (lfd);
          free (existing_lockfile);
          continue;
        }
      }
      else
        fail ("Process %d on host %s was dispatching this experiment.\n"
              "If that process no longer exists, you should delete\n"
              "%s (as seen from %s)\n"
              "and then you should delete\n"
              "%s (as seen from %s)\n"
              "if it contains the string\n"
              "%s\n"
              "(grep %s %s &> /dev/null && rm %s) || echo failed\n",
              lf_pid, lf_hostname, existing_lockfile, lf_hostname, lockfile_name, hostname (), existing_lockfile,
              existing_lockfile, lockfile_name, lockfile_name);
    }
    else if (s.st_nlink != 2)
      fail ("%s has %lu links!  That should not be possible.\n", filename, (long unsigned) s.st_nlink);
    else                        /* s.st_nlink == 2 */
      return 0;
  }
  return 0;
}
