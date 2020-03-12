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

//g++  -Wall -std=gnu++11 -I. -g -O2 -o integrate integrate.cc
/* integrate.c */


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <config.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <cerrno>
#include <error.h>
using std::vector;
using std::ifstream;
using std::string;
using std::istringstream;
using std::cout;

#define TREALLOC(buf, n) buf = (char *) realloc (buf, (n) * sizeof *(buf))
#define DIE (error_at_line (1, 0, __FILE__, __LINE__, "fatal error"), false)

#define SF 25000
#define IEDTHZ 200
#define BETWEEN (SF / IEDTHZ)
#define IEDTINCRE (10000 / IEDTHZ)
#define DECAY (exp (-(1000.0 / SF) / 200))

#define BUFSZ 200000
int malloc_debug;

static short buf[BUFSZ];
static FILE *infile, *bin, *edt;
static int amin, amax, mean, action;
static double amean;
static bool copy0 = true;

int marker_id;
vector<int> marker;

enum {INTEGRATE, DIGITIZE};

#if defined i386 && defined __GNUC__
static inline short int
swap_bytes (short int x)
{
	__asm__("xchgb %b0,%h0"		/* swap bytes		*/
		: "=q" (x)
		:  "0" (x));
	return x;
}
#else
static inline short int
swap_bytes (short val)
{
  unsigned short uval;
  uval = *(unsigned short *)&val;
  uval = ((uval >> 8) & 0xff) | ((uval & 0xff) << 8);
  return *(short *)&uval;
}
#endif

#ifdef WORDS_BIGENDIAN
# define l_e_to_host(x) swap_bytes (x)
# define host_to_b_e(x) (x)
#else
# define l_e_to_host(x) (x)
# define host_to_b_e(x) swap_bytes (x)
#endif

static char *
change_filetype_trim (const char *file_name, const char *suffix, int trim)
{
  char *dot, *p;
  static char *new_file;

  TREALLOC (new_file, strlen (file_name) + strlen (suffix) + 1);
  strcpy (new_file, file_name);
  if ((dot = strrchr (new_file, '.')) == 0)
    strcat (new_file, suffix);
  else
    strcpy (dot, suffix);
  if (trim
      && ((p = strrchr (new_file, '/')) != 0
	  || (p = strrchr (new_file, '\\')) != 0))
    return p + 1;
  else
    return new_file;
}

static void
get_maxmin (void)
{
  int cnt, n, analog;
  double decay = DECAY;

  analog = 0;
  amin = SHRT_MAX;
  amax = SHRT_MIN;

  int32_t t = 0;
  bool copy = copy0;
  size_t marker_idx = 0;
  double sum = 0;
  int count = 0;
  while ((cnt = fread (buf, sizeof *buf, BUFSZ, infile)) > 0)
    for (n = 0; n < cnt; n++) {
      if (marker_idx < marker.size() && t >= marker[marker_idx]) {
        copy = !copy;
        marker_idx++;
      }
      ++t;
      if (copy) {
        sum += l_e_to_host (buf[n]);
        ++count;
      }
    }
  mean = sum / count;
  rewind (infile);

  t = 0;
  copy = copy0;
  marker_idx = 0;
  while ((cnt = fread (buf, sizeof *buf, BUFSZ, infile)) > 0)
    for (n = 0; n < cnt; n++) {
      if (marker_idx < marker.size() && t >= marker[marker_idx]) {
        copy = !copy;
        marker_idx++;
      }
      ++t;
      if (copy) {
        analog = (action == DIGITIZE
                  ? l_e_to_host (buf[n])
                  : analog * decay + abs (l_e_to_host (buf[n]) - mean));
        if (analog < amin) amin = analog;
        if (analog > amax) amax = analog;
        amean += analog;
      }
    }
  amean /= count;
  (feof (infile) && !ferror (infile)) || DIE;
  clearerr (infile);
}

static void
do_edt (int code)
{
  short val;
  int itim, iedtcntr, ianout, offset, analog = amean, cnt, n;
  double factor, decay = DECAY;

  factor = 4090.0 / (amax - amin);
  offset = -2047 - factor * amin;

  iedtcntr = 1;
  analog = amean;
  itim = 0;
  fprintf (edt, "%5d%10d\n", 33, 3333333);
  fprintf (edt, "%5d%10d\n", 33, 3333333);
  int32_t t = 0;
  bool copy = copy0;
  size_t marker_idx = 0;
  while ((cnt = fread (buf, sizeof *buf, BUFSZ, infile)) > 0) {
    for (n = 0; n < cnt; n++) {
      if (marker_idx < marker.size() && t >= marker[marker_idx]) {
        copy = !copy;
        marker_idx++;
      }
      ++t;
      
      buf[n] = host_to_b_e (val = l_e_to_host (buf[n]));
      if (copy) {
        analog = (action == DIGITIZE
                  ? val
                  : analog * decay + abs (val - mean));
      }
      else analog = amean;
      if (iedtcntr == BETWEEN) {
        if (copy) {
          if ((ianout = analog * factor + offset) < 0)
            ianout += 4096;
          ianout += code * 4096;
          fprintf (edt, "%5d%10d\n", ianout, itim);
        }
	iedtcntr = 0;
	itim = itim + IEDTINCRE;
      }
      iedtcntr++;
    }
    if (action == INTEGRATE)
      (int)fwrite (buf, sizeof *buf, cnt, bin) == cnt || DIE;
  }
}

static vector<int>
read_times (const char *filename, int marker_id)
{
  ifstream f(filename);
  if (!f)
    error (1, errno, "error opening %s", filename);
  vector<int> t;
  string line;
  while (getline (f, line)) {
    uint64_t time;
    int id;
    time = atoi (line.substr (5).c_str());
    istringstream (line.substr (0, 5)) >> id;
    if (id == marker_id)
      t.push_back (time * 25000 / 10000);
  }
  return t;
}

int
main (int argc, char **argv)
{
  int code;
  const char *i, *d, *istr = "integrate", *dstr = "digitize";

  i = strstr (argv[0], istr);
  d = strstr (argv[0], dstr);
  if (i && !d)
    action = INTEGRATE;
  else if (d && !i)
    action = DIGITIZE;
  else {
    printf ("program must have \"%s\" or \"%s\" in the name, but not both\n", istr, dstr);
    exit (1);
  }

  if (argc == 6) {
    if (string (argv[5]) == "flip")
      copy0 = false;
    --argc;
  }
  if (argc == 5) {
    marker_id = atoi (argv[3]);
    marker = read_times (argv[4], marker_id);
  }
  else if (argc != 3) {
    cout << R"(
  usage: integrate FILE.chan ANALOG_ID [CUTCODE CUT_EDT [flip]]

  Rectifies and integrates the signal in FILE.chan and writes the
  result to FILE.edt with the specified ANALOG_ID.  Also writes
  FILE.bin, which is a big-endian copy of FILE.chan.  If
  FILE.edt or FILE.bin exists, it will be overwritten.  Also
  appends an I to FILE.status.
  
  If CUTCODE and CUT_EDT are specified, regions will be left out of
  the integrated signal in FILE.edt.  The regions to be left out are
  specified by an event code in an .edt file.  CUTCODE is the id of
  the event code, and CUT_EDT is the name of the .edt file.  The
  region between the first and second event of the given id is left
  out, and between the third and fourth, etc.  In general, the region
  between each odd numbered event (start marker) and the following
  even numbered event (stop marker) is left out.  If 'flip' is
  specifed, the regions that would have been left out are kept, and
  vice-versa. The content of the left-out regions has no effect on the
  integrated signal.

)";
    exit (1);
  }

  (infile = fopen (argv[1], "rb")) || DIE;
  (edt = fopen (change_filetype_trim (argv[1], ".edt", 0), "wb")) || DIE;
  if (action == INTEGRATE)
    (bin = fopen (change_filetype_trim (argv[1], ".bin", 0), "wb")) || DIE;
  code = atoi (argv[2]);

  get_maxmin ();

  rewind (infile);

  do_edt (code);

  {
    FILE *f;
    (f = fopen (change_filetype_trim (argv[1], ".status", 0), "a")) || DIE;
    putc (action == INTEGRATE ? 'I' : 'D', f);
    fclose (f);
  }

  return 0;
}
