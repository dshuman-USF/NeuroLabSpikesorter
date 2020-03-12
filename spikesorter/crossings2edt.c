#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
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
  fprintf (stderr, "usage: %s whatever.chan threshold unit_code\n", progname);
  exit (0);
}

int
main (int argc, char **argv)
{
  short val;
  int armed = 0;
  FILE *chan, *edt;
  int sample_count, code;
  int sample = 0;
  int lasttime = 0, now;

  progname = argv[0];
  if (argc < 4)
    usage ();

  (chan = fopen (argv[1], "rb")) || DIE;
  threshold = atoi (argv[2]);
  code = atoi (argv[3]);
  sample_count = filesize (chan) / sizeof (short);

  char *edtname = strdup (argv[1]);
  char *dot = strrchr (edtname, '.');
  if (!dot) usage ();
  if (strcmp (dot, ".chan") != 0) usage ();
  strcpy (dot, ".edt");
  (edt = fopen (edtname, "wb")) || DIE;

  fprintf (edt, "   33   3333333\n");
  fprintf (edt, "   33   3333333\n");

  while (fread (&val, sizeof val, 1, chan) == 1) {
    if ((now = time (0)) > lasttime) {
      lasttime = now;
      fprintf (stderr, "  %.0f%%\r", floor ((double)sample / sample_count * 100 + .5));
      fflush (0);
    }
    if (val < threshold)
      armed = 1;
    else if (val > threshold && armed) {
      armed = 0;
      fprintf (edt, "%5d%10d\n", code, (int)floor (sample / 2.5 + .5));
    }
    sample++;
  }
  return 0;
}

