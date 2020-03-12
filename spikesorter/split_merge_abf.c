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

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define DIE (fprintf (stderr, "error at line %d\n", __LINE__), perror (""), exit(1), 0)

#define BUFSIZE (1<<23)		/* must be >= 1024 and a multiple of 8 */
short inbuf[BUFSIZE], outbuf[BUFSIZE/16];
int start[] = {0,0,8,4,12,1,9,5,13,2,10,6,14,3,11,7,15};
double average;
FILE *of[17];
char *abf_file;

int strcasecmp( const char *s, const char *d)
{
	for (;;) {
		if (*s != *d) {
			if (tolower(*s) != tolower(*d))
				return *s - *d;
		} else if (*s == '\0')
			break;
		s++;
		d++;
	}
	return 0;
}

void process (FILE *ifile)
{
	int readcount, from, to;
	int channel;
	
	fseek (ifile, 2048, SEEK_SET) && DIE;				//seek past the 2K header
	while ((readcount = fread (inbuf, sizeof *inbuf, BUFSIZE, ifile)) > 0) {
		(readcount % 16 == 0) || DIE;
		for (channel = 1; channel <= 16; channel++){
			if(of[channel] == 0)
				continue;
			for (from = start[channel], to = 0; from < readcount; from += 16, to++)
				outbuf[to] = inbuf[from];
			fwrite (outbuf, sizeof *outbuf, readcount/16, of[channel]);
		}
	}
}

void open_files (int ccnt, char **argv)
{
	int nmlen, chan, n;
	char *nmbuf, *ext, *p;

	if ((p = strrchr (abf_file, '/')) == 0)
		if ((p = strrchr (abf_file, '\\')) == 0)
			p = abf_file;
	if (*p == '/' || *p == '\\')
		p++;
	nmlen = strlen (p);
	(nmbuf = malloc (nmlen+4)) || DIE;
	strcpy (nmbuf, p);
	ext = nmbuf + nmlen - 4;
	
	for (n = 0; n < ccnt; n++){
		chan = atoi(argv[n]);
		(chan > 0 && chan <= 16) || DIE;
		sprintf (ext, "%02d.chan",chan);
		(of[chan] = fopen (nmbuf, "wb")) || DIE;
	}
}

int main (int argc, char **argv)
{
	int nmlen;
	char *nmbuf, *ext, *oldext, *nmbuf2, *ext2;
	char seq = 'A';		//what about lower case??
	FILE *f;
	char *argv16[16] = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"};
	
	if (argc < 2) {
		printf ("usage: %s whatever.abf [channel]...\n", argv[0]);
		exit (0);
	}
	nmlen = strlen (argv[1]);
	nmbuf = malloc (nmlen+2);
	nmbuf2 = malloc (nmlen+2);
	strcpy (nmbuf, argv[1]);
	strcpy (nmbuf2, argv[1]);
	ext = nmbuf + nmlen - 4;
	ext2 = nmbuf2 + nmlen - 4;
	oldext = strdup (ext);
	
	if (nmlen < 5 || strcasecmp (ext, ".abf") != 0) {
		printf ("filename must end in \".abf\"\n");
		exit (1);
	}
	
	abf_file = argv[1];
	
	if (argc == 2)		//open all 16 output files
		open_files (16, argv16);
	else						//open select files (<16)
		open_files (argc-2, argv+2);
	
	while ((f = fopen (nmbuf, "rb")) || (f = fopen (nmbuf2, "rb"))) {
		printf ("%s\n", nmbuf);
		process (f);
		fclose (f);
		sprintf (ext, "%c%s", seq, oldext);
		sprintf (ext2, "%c%s", tolower(seq++), oldext);
	}
	
	return 0;
}
// -*- c-basic-offset: 4; tab-width: 4; -*-
