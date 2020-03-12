#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "util.h"
typedef struct {double loc; int id; int uidx;} Entry;

char *progname;

void
usage (void)
{
  fprintf (stderr, "usage: %s xxx.edt > yyy.pos\n", progname);
  exit (0);
}

int
get_data (FILE *edt, int *code, double *time)
{
  static char *line;
  static size_t len;
  int chars_read;

  if ((chars_read = getline (&line, &len, edt)) < 0)
    return 0;
  chars_read > 15 || DIE;
  sscanf (line + 5, "%lf", time) == 1 || DIE;
  (*time >= 0 && *time < 9999999999.0 && floor (*time) == *time) || DIE;

  line[5] = 0;
  sscanf (line, "%d", code) == 1 || DIE;
  (*code > 0 && *code < 99999) || DIE;
  return 1;
}

int
main (int argc, char **argv)
{
  FILE *edt;
  Entry etmp;

  progname = argv[0];
  if (argc != 2)
    usage ();

  (edt = fopen (argv[1], "rb")) || DIE;

  int code;
  double time;
  get_data (edt, &code, &time);
  get_data (edt, &code, &time);
  while (get_data (edt, &code, &time)) {
    etmp.loc = time * 2.5;
    etmp.id = code;
    etmp.uidx = 0;
    fwrite (&etmp, sizeof etmp, 1, stdout) == 1 || DIE;
  }
  return 0;
}

