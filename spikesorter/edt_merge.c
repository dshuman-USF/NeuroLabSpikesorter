/* Copyright 2005-2020 Kendall F. Morris

    This file is part of the Spiksorter software suite.

    The Spikesorter software suite is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation, either
    version 3 of the License, or (at your option) any later version.

    The suite is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the suite.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <config.h>
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <search.h>
#include "getline.h"

#ifndef HAVE_ASPRINTF
int asprintf(char **buffer, char *fmt, ...);
#endif

#define DIE (fprintf (stderr, "fatal error in %s at line %d \n", __FILE__, __LINE__), exit (1), 0)
#define TMALLOC(buf, n) (buf = malloc ((n) * sizeof *(buf))) || DIE
#define TCALLOC(buf, n) (buf = calloc ((n), sizeof *(buf))) || DIE
#define TREALLOC(buf, n) (buf = realloc (buf, (n) * sizeof *(buf))) || DIE
#define TMEMCPY(to, from, n) memcpy ((to), (from), n*sizeof *(to))
#define TMEMSET(to, from, n) memset ((to), (from), n*sizeof *(to))
#undef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#undef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))

typedef struct {
  int channel;
  int code_start;
  int code_end;
  int array_size;
  int array_electrode;
  int unit_count;
  int last_code;
  int next_code;
  double last_time, next_time;
  int *unit;
  int *num_map;
  FILE *edt;
} ChannelInfo;

#define MAX_DIG_CODE 999
 
int
first_choice (ChannelInfo *c, int unit_idx)
{
  int code;
  int rounded_array_size = (int)(ceil (c->array_size / 10) * 10);

  // unit_idx is zero-based
  // array_electrode is one-based

  code = (c->code_start + (unit_idx * rounded_array_size) + (c->array_electrode - 1));
  if (code > c->code_end)
    return 0;
  return code;
}

int
second_choice (ChannelInfo *c, char *taken)
{
  int n;
  for (n = c->code_end; n >= c->code_start; n--)
    if (!taken[n])
      return n;
  exit (DIE);
}

static char *base_name;
static int total_unit_count;

void
get_merge_nums (ChannelInfo *channel_info, int channel_count, int max_unit_count, char *taken)
{
  int unit_idx, chanidx, merge_num;
  ChannelInfo *c;
  char *file_name;
  char *file_name_2;
  FILE *f;
  FILE *f2;
  ENTRY e;

  if (asprintf (&file_name  , "%s.rpt", base_name) == -1) exit (1);
  if (asprintf (&file_name_2, "%s.nam", base_name) == -1) exit (1);
  (f  = fopen (file_name  , "w")) || DIE;
  (f2 = fopen (file_name_2, "w")) || DIE;

  hcreate (total_unit_count) || DIE;

  taken[0] = 1;
  for (unit_idx = 0; unit_idx < max_unit_count; unit_idx++) {
    for (chanidx = 0; chanidx < channel_count; chanidx++) {
      c = channel_info + chanidx;
      if (unit_idx >= c->unit_count)
	continue;
      (merge_num = first_choice (c, unit_idx)) <= MAX_DIG_CODE || DIE;
      if (taken[merge_num])
	(merge_num = second_choice (c, taken)) <= MAX_DIG_CODE || DIE;
      taken[merge_num] = 1;
      c->num_map[c->unit[unit_idx]] = merge_num;
      fprintf (f, "channel %2d (electrode %2d) unit %2d mapped to %3d\n",
	       c->channel, c->array_electrode, c->unit[unit_idx], merge_num);
      {
	char *name_string;
	char *name = "IRSVTPFDZA";
	char namc;
	char suffix[2];
	suffix[1] = 0;
	suffix[0] = unit_idx ? ('a' + unit_idx - 1) : 0;
	(c->code_start >= 1 && c->code_start < 1000) || DIE;
	namc = name[(c->code_start - 1) / 100];
	while (1) {
	  ENTRY *ep;
	  if (asprintf (&e.key, "%c%d%s", namc, c->array_electrode, suffix) == -1) exit (1);
	  if ((ep = hsearch (e, FIND)) == 0)
	    break;
	  free (e.key);
	  if (suffix[0] == 0) suffix[0] = 'a'; else suffix[0]++;
	  suffix[0] <= 'z' || DIE;
	}
	name_string = e.key;
	e.key = strdup (name_string);
	hsearch (e, ENTER) || DIE;
	fprintf (f2, "%d %s %d\n", c->channel, name_string, merge_num);
	free (name_string);
      }
    }
  }
  fclose (f);
  fclose (f2);
  free (file_name);
}

#define HDR_STR "   33   3333333\n"
void
get_header (ChannelInfo *c)
{
  static char *line;
  static size_t len;
  int n, i;
  
  for (i = 0; i < 2; i++) {
    if (getline (&line, &len, c->edt));
    strncmp (line, HDR_STR, 15) == 0 || DIE;
    n = strlen (line);
    n == 16 || n == 17 || DIE;
  }
}

int
get_data (ChannelInfo *c, int *code, double *time)
{
  static char *line;
  static size_t len;
  int chars_read;

  if ((chars_read = getline (&line, &len, c->edt)) < 0)
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
min_time (ChannelInfo *channel_info, int channel_count, int *code, double *time)
{
  int n, min;
  ChannelInfo *c = channel_info;

  for (min = 0, n = 1; n < channel_count; n++)
    if (c[n].next_time < c[min].next_time)
      min = n;
  if (c[min].next_time == DBL_MAX)
    return 0;
  if ((n = c[min].next_code) <= c[min].unit[c[min].unit_count-1] && c[min].num_map[n])
    n = c[min].num_map[n];
  *code = n;
  *time = c[min].next_time;
  if (c[min].last_time > c[min].next_time) {
    printf ("Codes are out of order in input .edt file\n");
    printf ("channel %d, code %d time %.17g last %.17g\n", min, *code, *time, c[min].last_time);
    exit (DIE);
  }
  c[min].last_time = c[min].next_time;
  c[min].last_code = c[min].next_code;
  if (!get_data (c + min, &c[min].next_code, &c[min].next_time))
    c[min].next_time = DBL_MAX;
  return 1;
}

static int noedt;

static void
merge (void)
{
  int channel_count, n, i, max_unit_count, code;
  double time;
  ChannelInfo *channel_info, *c = 0;
  static char taken[MAX_DIG_CODE + 1];
  char *file_name;
  FILE *out;
  FILE *f;
  
  if (asprintf (&file_name, "%s.chs", base_name) == -1) exit (1);
  (f = fopen (file_name, "r")) || DIE;
  fscanf (f, "%d", &channel_count) == 1 || DIE;
  TMALLOC (channel_info, channel_count);
  max_unit_count = 0;
  for (n = 0; n < channel_count; n++) {
    c = channel_info + n;
    fscanf (f, "%d %d %d %d %d %d",
	    &c->channel, &c->code_start, &c->code_end, &c->array_size, &c->array_electrode, &c->unit_count) == 6 || DIE;
    total_unit_count += c->unit_count;
    if (c->unit_count > max_unit_count)
      max_unit_count = c->unit_count;
    TMALLOC (c->unit, c->unit_count);
    for (i = 0; i < c->unit_count; i++)
      fscanf (f, "%d", &c->unit[i]) == 1 || DIE;
    c->unit[c->unit_count-1] <= MAX_DIG_CODE || DIE;
    TCALLOC (c->num_map, c->unit[c->unit_count-1] + 1);
  }
  get_merge_nums (channel_info, channel_count, max_unit_count, taken);
  if (noedt)
    exit (0);
  for (n = 0; n < channel_count; n++) {
    c = channel_info + n;
    if (asprintf (&file_name, "%s_%02d.edt", base_name, c->channel) == -1) exit (1);
    (c->edt = fopen (file_name, "r")) || DIE;
    free (file_name);
    get_header (c);
    c->last_code = 0;
    c->last_time = 0;
    if (!get_data (c, &c->next_code, &c->next_time))
      c->next_time = DBL_MAX;
  }
  if (asprintf (&file_name, "%s.edt", base_name) == -1) exit (1);
  (out = fopen (file_name, "w")) || DIE;
  free (file_name);
  fputs (HDR_STR, out);
  fputs (HDR_STR, out);
  while (min_time (channel_info, channel_count, &code, &time))
    fprintf (out, "%5d%10.0f\n", code, time);
}

int
main (int argc, char **argv)
{
  if (argc == 3 && strcmp (argv[2], "noedt") == 0)
    noedt = 1;
  else
    argc == 2 || DIE;
  base_name = argv[1];
  merge ();
  return 0;
}
