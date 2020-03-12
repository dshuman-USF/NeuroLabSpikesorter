#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <errno.h>
#include <error.h>
#include <stdbool.h>
#include <string.h>
#include "util.h"
typedef struct {double loc; int id; int uidx;} Entry;

static long
filesize (FILE *f)
{
  long bytes;

  fseek (f, 0, SEEK_END);
  bytes = ftell (f);
  rewind (f);
  return bytes;
}

int threshold = -16384;
char *progname;

void
usage (void)
{
  fprintf (stderr, "usage: %s [-m THR2] POS CHAN THR CLUSTER... > OUTPOS\n", progname);
  exit (0);
}

bool
thr2ok (short *buf, int i, int n, int thr1, int thr2)
{
  if (thr2 >= 32767)
    return true;
  for ( ; i < n; ++i) {
    if (buf[i] > thr2)
      return false;
    if (buf[i] < thr1)
      return true;
  }
  return false;
}

int
get_thr2 (int *argc, char **argv)
{
  if (strcmp (argv[1], "-m") != 0)
    return 32768;
  int thr2 = atoi (argv[2]);
  for (int i = 3; i <= *argc; ++i)
    argv[i - 2] = argv[i];
  *argc -= 2;
  return thr2;
}

int
main (int argc, char **argv)
{
  short val;
  int armed = 0;
  FILE *pos;
  FILE *chan;
  int count, sample_count, max, n;
  Entry *spike;
  int spikeidx = 0, sample = 0, id0;
  int lasttime = 0, now;
  int *drop;
  Entry etmp;

  progname = argv[0];

  if (argc < 5)
    usage ();
  int thr2 = get_thr2 (&argc, argv);
  if (argc < 5)
    usage ();

  (pos = fopen (argv[1], "rb")) || DIE;
  (chan = fopen (argv[2], "rb")) || DIE;
  threshold = atoi (argv[3]);
  max = 0;
  for (n = 4; n < argc; n++) {
    if (atoi (argv[n]) < 1)
      usage ();
    if (atoi (argv[n]) > max)
      max = atoi (argv[n]);
  }
  TCALLOC (drop, max + 1);
  for (n = 4; n < argc; n++)
    drop[atoi (argv[n])] = 1;
  id0 = atoi (argv[4]);
  count = filesize (pos) / sizeof (Entry);
  sample_count = filesize (chan) / sizeof (short);
  filesize (pos) == count *  sizeof (Entry) || DIE;
  TMALLOC (spike, count);
  fread (spike, sizeof *spike, count, pos) == count || DIE;

  short *buf = malloc (sample_count * sizeof (short));
  if (fread (buf, sizeof val, sample_count, chan) != sample_count)
    error (1, errno, "error reading %s", argv[2]);
  for (int i = 0; i < sample_count; ++i) {
    val = buf[i];
    if ((now = time (0)) > lasttime) {
      lasttime = now;
      fprintf (stderr, "  %.0f%%\r", floor ((double)sample / sample_count * 100 + .5));
      fflush (0);
    }
    if (val < threshold)
      armed = 1;
    else if (val > threshold && armed) {
      armed = 0;
      if (thr2ok (buf, i, sample_count, threshold, thr2)) {
        for (; spikeidx < count && spike[spikeidx].loc < sample; spikeidx++)
          if (spike[spikeidx].id > max || !drop[spike[spikeidx].id])
            fwrite (&spike[spikeidx], sizeof *spike, 1, stdout) == 1 || DIE;
        etmp.loc = sample;
        etmp.id = id0;
        etmp.uidx = 0;
        fwrite (&etmp, sizeof etmp, 1, stdout) == 1 || DIE;
      }
    }
    sample++;
  }
  for (; spikeidx < count; spikeidx++)
    if (spike[spikeidx].id > max || !drop[spike[spikeidx].id])
      fwrite (&spike[spikeidx], sizeof *spike, 1, stdout) == 1 || DIE;
  return 0;
}

