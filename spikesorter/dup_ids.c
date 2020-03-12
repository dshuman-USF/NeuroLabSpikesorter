#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <error.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

typedef struct
{
  int code;
  int time;
} Spike;

static int map[1000][3];
static int id1, id2;
static char *fmt, *hdr;

int
main (int argc, char **argv)
{
  if (argc != 3) {
    printf (
            "\n"
            "usage: %s map_file edt_file > out_edt_file\n"
            "\n"
            "changes each id specified in map_file that it finds in edt_file\n"
            "to the two id's specified in the first line of map_file,\n"
            "with the offsets specified for that id.\n"
            "The first line should contain two id's with a space between.\n"
            "Each of the other lines should contain an id, and two offsets in ticks\n"
            "separated by spaces.\n"
            "The result is sent to stdout.\n"
            "\n"
            , argv[0]);
    return 0;
  }

  FILE *mapfile = fopen (argv[1], "r");
  if (!mapfile) error (1, errno, "can't open %s for read", argv[1]);
  FILE *edt = fopen (argv[2], "r");
  if (!edt) error (1, errno, "can't open %s for read", argv[2]);
  if (fscanf (mapfile, " %d %d", &id1, &id2) != 2)
    error (1, errno, "Can't read first line of map_file");
  if (id1 < 0 || id1 > 99999 || id2 < 0 || id2 > 99999)
    error (1, 0, "ids in first line of map_file out of range");
  int code, offset1, offset2, count = 0;
  while (fscanf (mapfile, " %d %d %d", &code, &offset1, &offset2) == 3) {
    if (code < 0 || code > 999)
      error (1, 0, "id out of range in map file");
    map[code][0] = 1;
    map[code][1] = offset1;
    map[code][2] = offset2;
    count++;
  }
  fprintf (stderr, "%d ids to duplicate\n", count);

  char line[17];
  if (fgets (line, 17, edt));
  if (fgets (line, 17, edt));
  hdr = strdup (line);
  int len = 0;
  if (line[15] == '\n')
    len = 16, fmt = "%5d%10d\n";
  else if (line[13] == '\n')
    len = 14, fmt = "%5d%8d\n";
  else
    error (1, 0, "Bad line length for edt file");

  struct stat s;
  stat (argv[2], &s);
  if (s.st_size % len != 0)
    error (1, 0, "edt file size is not a multiple of the line length");
  int in_count = s.st_size / len - 2;
  fprintf (stderr, "%d events in edt file\n", in_count);
  Spike *out = malloc (2 * in_count * sizeof *out);

  int o = 0;
  for (int i = 0; i < in_count; i++) {
    if (!fgets (line, 17, edt))
      error (1, errno, "error reading line %d of the edt file", i + 3);
    int time = atoi (line + 5);
    line[5] = 0;
    int code = atoi (line);
    if (code == 0)
      error (1, 0, "zero code");
    if (code >= 0 && code < 1000 && map[code][0]) {
      out[o].code = id1;
      out[o].time = time + map[code][1];
      o++;
      out[o].code = id2;
      out[o].time = time + map[code][2];
      o++;
    }
    else {
      out[o].code = code;
      out[o].time = time;
      o++;
    }
  }
  
  int out_count = o;


  int
  compare_spikes (const void *p1, const void *p2)
  {
    const Spike *s1 = (const Spike *)p1;
    const Spike *s2 = (const Spike *)p2;
    if (s1->time == s2->time)
      return (s1->code > s2->code) - (s1->code < s2->code);
    return (s1->time > s2->time) - (s1->time < s2->time);
  }
  
  qsort (out, out_count, sizeof out[0], compare_spikes);
  
  printf ("%s", hdr);
  printf ("%s", hdr);
  for (int i = 0; i < out_count; i++)
    printf (fmt, out[i].code, out[i].time);
  return 0;
}

// /raid/datamax/2004-07-29/clean/2004-07-29-ch25-corrections.txt

// ls /raid/sarahtemp/2004-07-29_001/2004-07-29_001_v61.edt

