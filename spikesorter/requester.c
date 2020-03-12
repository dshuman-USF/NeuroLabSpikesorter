#define _GNU_SOURCE
#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <dirent.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <limits.h>
#include <ctype.h>
#include <sys/file.h>
#include <sys/wait.h>
#include <sys/ptrace.h>
#include <linux/ptrace.h>
#include <sys/utsname.h>

static int
err (char *file, int line, int p, char *fmt, ...)
{
  FILE *f;
  va_list ap;
  char *msg;

  f = fopen ("/tmp/err", "a");
  fprintf (f, "%s line %d: ", file, line);
  va_start (ap, fmt);
  if (vasprintf (&msg, fmt, ap) == -1) exit (1);
  va_end (ap);
  if (p)
    fprintf (f, "%s: %s\n", msg, strerror (errno));
  else
    fprintf (f, "%s\n", msg);
  fclose (f);
  exit (1);
  return 1;
}
#define  ERR(fmt...) err (__FILE__, __LINE__, 0, ## fmt)
#define PERR(fmt...) err (__FILE__, __LINE__, 1, ## fmt)

static int
die (char *file, int line, int p, char *fmt, ...)
{
  va_list ap;
  char *msg;

  fprintf (stderr, "%s line %d: ", file, line);
  va_start (ap, fmt);
  if (vasprintf (&msg, fmt, ap) == -1) exit (1);
  va_end (ap);
  if (p)
    perror (msg);
  else
    fprintf (stderr, "%s\n", msg);
  exit (1);
  return 1;
}
#define  DIE(fmt...) die (__FILE__, __LINE__, 0, ## fmt)
#define PDIE(fmt...) die (__FILE__, __LINE__, 1, ## fmt)

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
  if (getline (&line, &len, statfile) == -1) PDIE ("");
  p = strrchr (line, ')') + 2;
  if (sscanf (p, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %lu", &starttime) != 1)
    return 0;                   /* process must just have died */
  fclose (statfile);
  free (statfilename);
  return starttime;
}

static bool tail_done;

static void
set_tail_done (int signum)
{
  tail_done = 1;
}

static void
ignore (int signum)
{
}

static char *nohup_filename;
static char *job_filename;

static void
tail (void)
{
  FILE *f;
  static char *line = NULL;
  static size_t len = 0;
  bool done = false;
  sigset_t mask, oldmask;

  signal(SIGUSR1, set_tail_done);
  signal (SIGALRM, ignore);
  if ((f = fopen (nohup_filename, "r")) == NULL)
    PDIE ("Can't open %s", nohup_filename);
  fseek (f, 0, SEEK_END);
  if (ftell (f) < 256)
    fseek (f, 0, SEEK_SET);
  signal (SIGALRM, ignore);
  sigemptyset (&mask);
  sigaddset (&mask, SIGUSR1);
  sigaddset (&mask, SIGALRM);

  while (!done) {
    if (tail_done)
      done = true;
    while (getline (&line, &len, f) > 0)
      printf ("%s", line);
    clearerr (f);
    sigprocmask (SIG_BLOCK, &mask, &oldmask);
    if (!tail_done) {
      alarm (1);
      sigsuspend (&oldmask);
    }
    sigprocmask (SIG_UNBLOCK, &mask, NULL);
  }
  fclose (f);
  exit (0);
}

static void
run_job (void)
{
  static char *line = NULL;
  static size_t len = 0;
  pid_t jpid, tpid;
  int pfd[2];

  printf ("requester: ready\n");

  if (getline (&line, &len, stdin) < 1)
    return;

  if (pipe(pfd) == -1) PDIE ("pipe");
  jpid = fork ();
  if (jpid) {
    char buf;
    close (pfd[1]);
    if (read (pfd[0], &buf, 1));     /* wait for nohup.out to be created */
    close (pfd[0]);
    if ((tpid = fork ())) {
      int status;
      waitpid (jpid, &status, 0);
      kill (tpid, SIGUSR1);
      waitpid (tpid, NULL, 0);
    }
    else tail ();
  }
  else {
    int fd;
    char *cmd;

    close (pfd[0]);

    if (setsid () == -1) PDIE ("");
    if (close (0) !=  0) PERR ("");
    if (close (1) !=  0) PERR ("");
    if (close (2) !=  0) PERR ("");

    if ((fd = open ("/dev/null", O_RDONLY)) != 0) PERR ("%d", fd);
    if ((fd = open (nohup_filename, O_WRONLY|O_CREAT|O_TRUNC, 0744)) != 1) PERR ("%d", fd);
    
    if ((fd = dup (1)) != 2) PERR ("%d", fd);
    close (pfd[1]);
    if (asprintf (&cmd, "{ %s\n}", line) == -1) exit (1);
    if (execl ("/bin/bash", "/bin/bash", "-c", cmd, NULL) == -1) PDIE ("");
  }
}

static int
check_exe (pid_t pid, char *path)
{
  char *exe_name;
  int linklen, needed;
  char *delstr = " (deleted)";
  int dellen = strlen (delstr);
  static char *buf;
  static int bufsize;
  
  needed = strlen (path) + dellen + 1;
  if (bufsize < needed)
    buf = realloc (buf, bufsize = needed);

  if (asprintf (&exe_name, "/proc/%d/exe", pid) == -1) exit (1);
  linklen = readlink (exe_name, buf, bufsize);
  free (exe_name);
  if (linklen < 0)
    return 0;
  buf[linklen] = 0;
  if (linklen >= dellen && strcmp (buf + linklen - dellen, delstr) == 0)
    buf[linklen - dellen] = 0;

  //  printf ("pid %d exe is %s\n", pid, buf); /* debug */
  
  return strcmp (path, buf) == 0;
}

static int
has_open (pid_t pid, void *args)
{
  char **file = (char **) args;
  char *dirname;
  DIR *dir;
  struct dirent *de;
  char *link_buf;
  int len, linklen, maxlen, count;
  int found = 0, n;

  if (file[0] == 0 || !check_exe (pid, file[0]))
    return 0;

  count = maxlen = 0;
  for (n = 0; file[n]; n++)
    if ((len = strlen (file[n])) > maxlen)
      maxlen = len;
  count = n - 1;
  link_buf = calloc (maxlen + 2, sizeof *link_buf);

  if (asprintf (&dirname, "/proc/%d/fd", pid) == -1) exit (1);
  
  if ((dir = opendir (dirname)) == NULL && errno != EACCES)
    PDIE ("Can't open %s", dirname);
  errno = 0;
  if (dir != NULL)
    while ((errno = 0, (de = readdir (dir)) != NULL))
      if (de->d_type == DT_LNK) {
        int fd;
        char *path;

        fd = atoi (de->d_name);
        if (fd < 0 || fd >= count)
          break;
        if (asprintf (&path, "%s/%s", dirname, de->d_name) == -1) exit (1);
        if ((linklen = readlink (path, link_buf, maxlen + 1)) == -1)
          PERR ("Can't read link %s\n", path);
        free (path);
        link_buf[linklen] = 0;
        if (linklen > maxlen || strcmp (link_buf, file[1 + fd]) != 0)
          break;
        found++;
      }
  if (de == NULL && errno != 0)
    PDIE ("Error reading %s", dirname);
  closedir (dir);
  free (dirname);
  free (link_buf);
  return found == count;
}

static pid_t
find_pid (int test (int pid, void *args), void *args)
{
  DIR *dir;
  struct dirent *de;
  char *procdir = "/proc";
  long unsigned starttime, earliest_starttime = ULONG_MAX;
  pid_t earliest_pid = 0;

  if ((dir = opendir (procdir)) == NULL)
    PDIE ("Can't open %s", procdir);
  while ((errno = 0, (de = readdir (dir)) != NULL))
    if (de->d_type == DT_DIR) {
      int n;
      pid_t pid;
      int len = strlen (de->d_name);
      for (n = 0; n < len; n++)
        if (!isdigit (de->d_name[n]))
          break;
      if (n == len) {
        pid = atoi (de->d_name);
        if (test (pid, args) && (starttime = get_starttime (pid)) < earliest_starttime) {
          earliest_starttime = starttime;
          earliest_pid = pid;
        }
      }
    }
  if (errno != 0)
    PDIE ("Error reading %s", procdir);
  return earliest_pid;
}

static char *lockfile_name;

static pid_t
job_running (void)
{
  char *nohup_path;
  char *file[] = {"/bin/bash", "/dev/null", "", "", lockfile_name, 0};
  enum {EXE, FD0, FD1, FD2, FCNT};
  pid_t pid;

  if (FCNT >= sizeof file / sizeof file[0]) DIE("");
  
  nohup_path = canonicalize_file_name (nohup_filename);
  if (nohup_path == NULL) {
    if (errno == ENOENT)
      return 0;
    else
      PDIE ("Can't canonicalize %s", nohup_filename);
  }

# define FILL(fd,f) if (strcmp (file[fd], "") != 0) DIE (""); file[fd] = f

  FILL (FD1, nohup_path);
  FILL (FD2, nohup_path);

  pid = find_pid (has_open, file);
  free (nohup_path);
  return pid;
}

static int
restart (pid_t pid, int status)
{
  siginfo_t siginfo;
  int sig = WSTOPSIG (status);

  if (ptrace (PTRACE_SETOPTIONS, pid, 0, PTRACE_O_TRACEEXEC) != 0)
    PDIE ("PTRACE_SETOPTIONS");
  if (ptrace (PTRACE_GETSIGINFO, pid, 0, &siginfo) == 0) {
    if (siginfo.si_code == ((PTRACE_EVENT_EXEC << 8) | SIGTRAP)
        || sig == SIGSTOP || sig == SIGTSTP || sig == SIGTTIN || sig == SIGTTOU)
      sig = 0;
    ptrace (PTRACE_CONT, pid, 0, sig);
  }
  return 1;
}

static off_t
filesize (FILE *f)
{
  struct stat s;
  
  if (fstat (fileno (f), &s) != 0) PDIE ("");
  return s.st_size;
}

static void
print_file (char *filename)
{
  FILE *f;
  int c;
  if ((f = fopen (filename, "r"))) {
    while ((c = getc (f)) != EOF)
      putchar (c);
    fclose (f);
  }
}

static bool
is_job (char *job)
{
  FILE *f;
  int job_number;
  if ((f = fopen (job_filename, "r")) == NULL)
    return false;
  if (fscanf (f, "%d", &job_number) != 1)
    DIE ("job_number file contents do not start with a number");
  fclose (f);
  return atoi (job) == job_number;
}

static void
print_status_file (void)
{
  FILE *f;
  char *buf;
  off_t sz;
  char *nm;

  if ((f = fopen (job_filename, "r")) == NULL)
    return;
  if ((buf = malloc (sz = filesize (f))) == NULL) PDIE ("malloc");
  if (fread (buf, 1, sz, f) != sz) PDIE ("fread");
  fclose (f);
  if (strtok (buf, " ") != NULL && (nm = strtok (NULL, "\r\n")) != NULL)
    print_file (nm);
  free (buf);
}

static void
print_done (void)
{
  printf ("requester: done ");
  print_status_file ();
  printf ("\n");
}

static void
wait_and_tail (pid_t jpid)
{
  pid_t tpid;

  if (ptrace (PTRACE_ATTACH, jpid, 0, 0) != 0) {
    if (errno == ESRCH)
      return;
    else
      PDIE ("PTRACE_ATTACH");
  }
  if ((tpid = fork ())) {
    int status;
    do {
      if (waitpid (jpid, &status, 0) == -1)
        PDIE ("waitpid");
    } while (WIFSTOPPED (status) && restart (jpid, status));

    kill (tpid, SIGUSR1);
    waitpid (tpid, NULL, 0);
    print_done ();
  }
  else tail ();
}

static char *
hostname (void)
{
  static struct utsname buf;
  if (uname (&buf) != 0) PDIE ("");
  return buf.nodename;
}

static void
print_status (char *job, char *filename, int cpu)
{
  printf ("requester: status ");
  if (job_running () && is_job (job))    
    printf ("running on %s_%d", hostname (), cpu);
  else 
    print_file (filename);
  printf ("\n");
}

int
main (int argc, char **argv)
{
  DIR *rqdir;
  FILE *rqfile;
  char *rqdirname;
  struct dirent *de;
  char *my_rq_filename;
  long long my_request_time;
  int fd;
  char *orig_dir;
  pid_t jpid;
  bool print_start_status = false;

  int old_requester (int argc, char **argv);
  
  if (argc == 2 || argc == 4)
    return old_requester (argc, argv);

  int cpu = atoi (argv[1]);

  if (asprintf (&rqdirname, "/var/local/requests/cpu%d", cpu) == -1) exit (1);
  if (asprintf (&lockfile_name, "%s/lock", rqdirname) == -1) exit (1);
  if (asprintf (&nohup_filename, "nohup_%s_%d.out", hostname (), cpu) == -1) exit (1);
  if (asprintf (&job_filename, "job_number_%s_%d", hostname (), cpu) == -1) exit (1);

  if (argc == 5 && strcmp (argv[2], "status") == 0) {
    print_status (argv[3], argv[4], cpu);
    return 0;
  }

  if (argc != 3) {
    printf ("usage: %s cpu experiment_name\n", argv[0]);
    return 0;
  }

  if (asprintf (&my_rq_filename, "%s/%s", rqdirname, argv[2]) == -1) exit (1);

  if ((orig_dir = get_current_dir_name()) == 0) PDIE("");
  rqfile = fopen (my_rq_filename, "r");
  if (rqfile) {
    int pid;
    long unsigned starttime;

    if (fscanf (rqfile, "%lld %d %lu %*s", &my_request_time, &pid, &starttime) != 3)
      DIE ("Can't parse request file %s", my_rq_filename);
    if (starttime == get_starttime (pid))
      if (kill (pid, SIGTERM) != 0)
        PDIE ("kill %d failed: %d\n", pid);
    fclose (rqfile);
  }
  else {
    struct timeval t;
    gettimeofday (&t, 0); 
    my_request_time = t.tv_sec * 1000000LL + t.tv_usec;
  }
  if ((rqfile = fopen (my_rq_filename, "w")) == NULL)
    DIE ("Can't open %s for write", my_rq_filename);
  fprintf (rqfile, "%16lld %d %lu %s\n", my_request_time, getpid (), get_starttime (getpid ()), argv[2]);
  fclose (rqfile);

  // no files open here, so lockfile gets fd == 3

  if ((fd = open (lockfile_name, O_RDONLY|O_CREAT, 0666)) < 0)
    PDIE ("Can't open %s", lockfile_name);
  if (fd != 3) DIE ("");

  printf ("requester: ");
  if ((jpid = job_running ()) != 0) {
    printf ("running ");
    print_file (job_filename);
    printf ("\n");
    wait_and_tail (jpid);
  }
  else {
    if (flock (fd, LOCK_EX|LOCK_NB) != 0)
      printf ("waiting\n");
    else
      print_start_status = true;
  }

  while (flock (fd, LOCK_EX) == 0) {
    int myturn  = 1;
    if (chdir (rqdirname) != 0)
      PDIE ("Can't change directory to %s", rqdirname);
    if ((rqdir = opendir (rqdirname)) == NULL)
      PDIE ("Can't open %s", rqdirname);
    errno = 0;
    while (errno = 0, (de = readdir (rqdir)) != NULL) {
      struct stat s;
      if (stat (de->d_name, &s) != 0) PDIE ("stat");
      if (S_ISREG (s.st_mode)) {
        //        long long time0, time1;
        int pid;
        long unsigned starttime;
        long long request_time;
        if (strcmp (de->d_name, argv[2]) == 0 || strcmp (de->d_name, "lock") == 0)
          continue;
        //        rdtscll (time0);
        if ((rqfile = fopen (de->d_name, "r")) == 0)
            PDIE ("Can't open %s for read", de->d_name);
        if (fscanf (rqfile, "%lld %d %lu %*s", &request_time, &pid, &starttime) != 3)
          DIE ("Can't parse request file %s", my_rq_filename);
        fclose (rqfile);
        //        rdtscll (time1);
        //        printf ("%lld\n", time1 - time0);
        
        if (starttime != get_starttime (pid))
          continue;               /* process not running */
        if (request_time < my_request_time) {
          myturn = 0;
          errno = 0;
          break;
        }
      }
    }
    if (de == NULL && errno != 0)
      PDIE ("Error reading %s", rqdirname);
    closedir (rqdir);
    if (chdir (orig_dir) != 0)
      PDIE ("Can't change directory to %s", orig_dir);
    
    if (myturn) {
      run_job ();
      print_done ();
    }
    else if (print_start_status) {
      printf ("waiting\n");
      print_start_status = false;
    }

    flock (fd, LOCK_UN);
    if (!myturn)
      sleep (1);
  }
  return 0;
}
