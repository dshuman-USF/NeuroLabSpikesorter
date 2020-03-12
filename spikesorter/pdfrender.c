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

/* pdfrender.c */

#include "nde.h"
#include <string.h>

#if HAVE_SPAWNVP
#include <process.h>
#endif

#define TH 25

static unsigned lastcolor;
static int color_set;

void start_pdf_page (FILE *f)
{
  fprintf (f,
	   "0 612 translate\n"
	   ".48 -.48 scale\n"
	   "25.5 37.5 translate\n"
	   "/Helvetica findfont [%d 0 0 -%d 0 0] makefont setfont\n",
	   TH, TH);
  color_set = 0;
}

/*
  /d {moveto 0 .01 rlineto stroke} bind def
*/

static char *pdf_rootname;
static int pdf_file_count;

static void
write_file_header (FILE *f)
{
  fprintf (f, "%s",
	   "<< /PageSize [792 612] /Orientation 0 >> setpagedevice \n"
	   "/m {moveto} bind def\n"
	   "/l {lineto} bind def\n"
	   "/s {stroke} bind def\n"
	   "/d {moveto closepath stroke} bind def\n"
	   );
}

FILE *
create_pdf_file (char *filename)	
{
  FILE *f;
  int len;

  if (filename == 0)
    filename = "psrender.ps";
  pdf_rootname == 0 || DIE;
  pdf_rootname = strdup (filename);
  (len = strlen (pdf_rootname)) > 3 || DIE;
  strcmp (pdf_rootname + len - 3, ".ps") == 0 || DIE;
  pdf_rootname[len - 3] = 0;
  pdf_file_count == 0 || DIE;
  pdf_file_count = 1;
  (f = fopen (filename, "w")) || DIE;
  return f;
}

void
resize_pdf (/*@unused@*/FILE *winp, /*@unused@*/ int W, /*@unused@*/ int H)
{
}

void
setcolor (FILE *f, unsigned color)
{
  float red, green, blue;
  red = (double) ((color >> 24)& 0xff) / 255;
  green = (double) ((color >> 16) & 0xff) / 255;
  blue = (double) ((color >> 8) & 0xff) / 255;

  if (red == green && green == blue)
    fprintf (f, "%.3f setgray\n", red);
  else
    fprintf (f, "%.3f %.3f %.3f setrgbcolor\n", red, green, blue);
  color_set = 1;
  lastcolor = color;
}

void
pdfrender (SnapList *sl, FILE *f, int imageW, int imageH, int x0, int y0, int scale, unsigned color)
{
#define scaleW(i) ((int)(((i) * ((imageW-1.0)/63.0))))
#define scaleH(i) ((imageH-1) - (int)(sl->snap[n]->val[i]/scale*(imageH/2) + (imageH/2)))
#define W 1600
#define H 1200
	int n, i, ax, ay, bx, by;
	float red, green, blue;

	if (!color_set || color != lastcolor) {
	  red = (double) ((color >> 24)& 0xff) / 255;
	  green = (double) ((color >> 16) & 0xff) / 255;
	  blue = (double) ((color >> 8) & 0xff) / 255;

	  if (red == green && green == blue)
	    fprintf (f, "%.3f setgray\n", red);
	  else
	    fprintf (f, "%.3f %.3f %.3f setrgbcolor\n", red, green, blue);
	  color_set = 1;
	  lastcolor = color;
	}

	fprintf (f, "0 setlinewidth\n");
	if (imageW+x0 > W || imageH+y0 > H) {	
		printf ("The image may not fit in the window.  It may have been truncated.\n");
	}
	for (n = 0; n < sl -> count; n++) {
		for (i = 0; i < 63; i++) {
			ax = x0 + scaleW(i);
			bx = x0 + scaleW(i+1);
			ay = y0 + scaleH(i);
			by = y0 + scaleH(i+1);
			if(ax > W-1 || ax > x0 + imageW-1)
				ax = (x0+imageW-1 >= W-1) ? W-1 : x0 + imageW-1;
			if(bx > W-1 || bx > x0 + imageW-1)
				bx = (x0+imageW-1 >= W-1) ? W-1 : x0 + imageW-1;
			if(ay > H-1 || ay > y0 + imageH-1)
				ay = (y0+imageH-1 >= H-1) ? H-1 : y0 + imageH-1;
			if(ay < y0)
				ay = y0;
			if(by > H-1 || by > y0 + imageH-1)
				by = (y0+imageH-1 >= H-1) ? H-1 : y0 + imageH-1;
			if(by < y0)
				by = y0;
			if (i == 0)
			  fprintf (f, "%d %d m\n", ax, ay);
			else
			  fprintf (f, "%d %d l\n", bx, by);
		}
		fprintf (f, "s\n");
	}
	return;
}

void pdftext (FILE *f, int x, int y, char *txt)
{
  float red, green, blue;
  unsigned color = lastcolor;

  red = (double) ((color >> 24)& 0xff) / 255;
  green = (double) ((color >> 16) & 0xff) / 255;
  blue = (double) ((color >> 8) & 0xff) / 255;

  if (red == green && green == blue)
    fprintf (f,
	     "(%s) %d %d 2 copy 3 add 4 index stringwidth pop -%d 1 setgray rectfill m %.3f setgray show\n",
	     txt, x, y, TH, red);
  else
    fprintf (f,
	     "(%s) %d %d 2 copy 3 add 4 index stringwidth pop -%d 1 setgray rectfill m %.3f %.3f %.3f setrgbcolor show\n",
	     txt, x, y, TH, red, green, blue);

//dup true charpath pathbbox 2 index sub exch 3 index sub exch 1 setgray rectfill 3 1 roll m 0 setgray show
}

void pdflines (FILE *f, int imageW, int imageH, int x0, int y0, unsigned color)
{
	int i, ax, ay, bx, by;
	float red, green, blue;

	if (color != lastcolor) {
	  red = (double) ((color >> 24)& 0xff) / 255;
	  green = (double) ((color >> 16) & 0xff) / 255;
	  blue = (double) ((color >> 8) & 0xff) / 255;

	  if (red == green && green == blue)
	  fprintf (f, "%.3f setgray\n", red);
	  else
	  fprintf (f, "%.3f %.3f %.3f setrgbcolor\n", red, green, blue);
	  color_set = 1;
	  lastcolor = color;
	}

	fprintf (f, "0 setlinewidth\n");
	if (imageW+x0 > W || imageH+y0 > H) {	
		printf ("The image may not fit in the window.  It may have been truncated.\n");
	}
	for (i = 0; i < 63; i++) {
	  ax = x0 + scaleW(i);
	  bx = x0 + scaleW(i);
	  ay = y0;
	  by = (y0+imageH-1 >= H-1) ? H-1 : y0 + imageH-1;
	  fprintf (f, "%d %d m %d %d l s\n", ax, ay, bx, by);
	}

	return;
}

static void
newfile (FILE *f)
{
  char *newname;

  pdf_file_count++;
  pdf_rootname || DIE;
  if (asprintf (&newname, "%s_%d.ps", pdf_rootname, pdf_file_count) == -1) exit (1);
  freopen (newname, "w", f) == f || DIE;
  free (newname);
  write_file_header (f);
}

void end_pdf_page (FILE *f)
{
  fprintf (f, "showpage\n");
  if (ftell (f) > 1024 MEG)
    newfile (f);
}

#if HAVE_SPAWNVP
typedef  /*@null@*/ const char * charstar;
void
ps2pdf (void)
{
  static /*@observer@*/ char *gs = "gs.exe";
  static charstar const argv []
	      = {"gs", "-dSAFER", "-dCompatibilityLevel=1.2", "-q", "-dNOPAUSE",
		 "-dBATCH", "-sDEVICE=pdfwrite", "-sOutputFile=psrender.pdf", "-c",
		 ".setpdfwrite", "-f", "psrender.ps", 0};

  spawnvp (_P_WAIT, gs, argv);
}

void
rm_ps (void)
{
  static /*@observer@*/ char *cmd = "rm.exe";
  static charstar const argv[] = {"rm", "psrender.ps", 0};

  spawnvp (_P_WAIT, cmd, argv);
}
#endif

void pdfdisplay (FILE *f)
{
  fclose (f);
  pdf_rootname = 0;
  pdf_file_count = 0;
  //  system ("ps2pdf psrender.ps; acroread psrender.pdf");
//  ps2pdf ();
//  rm_ps ();
  //  system ("ps2pdf psrender.ps; #rm psrender.ps");
}

void create_ach (SnapList *sl, FILE *f, int x0, int y0, float bin_width, int axis_length, /*@unused@*/ unsigned color)
{
	int i, j, max = 0, ach[100], bin;
	float window, ach_scale;
	double time1, time2;
	
	window = bin_width*10000.;			// ???????? THIS IS NOT RIGHT ?????????????
	//	printf ("entering create_ach\n");
	//	printf ("sl -> count = %d\n", sl->count);
	//	printf ("window = %.1f\n",window);

	for (i = 0; i < 100; i++)		// initialize the ach array
		ach[i] = 0;
	for (i = 0; i < sl->count; i++) {		// for each snapshot in the list
		time1 = sl->snap[i]->sample;
		for (j = i+1; j < sl->count && (sl->snap[j]->sample)-time1 <= window; j++) {
			time2 = sl->snap[j]->sample;
			bin = (int)((time2-time1)/100);
			if (bin >= 0 && bin < 100)
			  ach[bin]++;
//			ach[(int)(((time2-time1)/25)/100)]++;
		}
	}
	
	for (i = 0; i < 100; i++) {			// find the value of the maximum bin in ach[]
		if (max < ach[i])
			max = ach[i];
	}
	ach_scale = (float)(axis_length)/max;
	//	printf ("max = %d; axis_length = %d; ach_scale = %f\n", max, axis_length, ach_scale);
/*	for (i = 0; i<100; i++) {
		printf ("ach[%d] = %d\n",i,ach[i]);
	}
	PAUSE;
*/
	for (i = 0; i < 100; i++) {
		ach[i] = ach[i] * ach_scale;
	}
/*	for (i = 0; i<100; i++) {
		printf ("ach[%d] = %d\n",i,ach[i]);
	}
	PAUSE;
*/
	
	fprintf (f, "%d %d m %d %d l %d %d l s\n", x0, y0-axis_length, x0, y0, x0+axis_length, y0);	//draw the axes
	for (i = 0; i < 100; i++) {
		fprintf (f, "%d %d m %d %d l s\n", x0+i+1, y0, x0+i+1, y0-ach[i]);
	}
	pdftext(f, x0+50, y0+50, "ACH");
	return;
}

//	WINDH=1.0/BINW			!WINDH = fraction of total histogram
//	LAG=(DELT*WINDH)+1	!determine in which bin the target event will be placed and...

void create_timeline (SnapList *sl, FILE *f, int length, int x0, int y0)
{
	int i;
	double lscale, total_interval, time;
	//	printf ("start create_timeline\n");
	fprintf (f, "%d %d m %d %d l s\n", x0, y0, x0+length, y0);
	total_interval = sl->snap[(sl->count-1)]->sample;
	lscale = length / total_interval;
	//	printf ("total interval = %f; length = %d; lscale = %20.15f\n", total_interval, length, lscale);
	for (i = 0; i < sl->count; i++) {
		time = sl->snap[i]->sample;
		fprintf (f, "%d %d m %d %d l s\n", (int)(x0+time*lscale), y0+10, (int)(x0+time*lscale), y0-10);		// make the tic mark
	}
	return;
}

void create_variance_histograms (SnapList *sl, SnapList *noise, FILE *f, int x0, int y0)
//void create_variance_histograms (SnapList *sl, FILE *f, int x0, int y0)
{
	int i, j, n, m, x, y, ii, point[SNAPLEN][100], point_noise[100], maxmax=0, bin;
	double value, scale, tmp=0, mean[SNAPLEN], std_dev[SNAPLEN], binwidth[SNAPLEN], lo_end=0;
	double mean_noise=0, std_dev_noise, binwidth_noise=0, lo_end_noise=0;
	long double sum=0;
//	char *label;
//	extern  double unzapped_sd;
	
	//	printf ("start create_var_histo\n");
	start_pdf_page(f);

	memset (point, 0, sizeof point);				// initializes the point array
	memset (point_noise, 0, sizeof point_noise);	// initializes the point_noise array
	memset (mean, 0, sizeof mean);					// initializes the mean array
//	to put labels on the histograms, create a string array (labels[SNAPLEN]), allocate it, load it with
//		sprintf commands
	//	printf ("fill array point\n");
	//	printf ("sl->count = %d\n", sl->count);
	//	printf ("sl->count = %d;  noise->count = %d\n", sl->count,noise->count);
//	PAUSE;
	//	printf ("sl->snap[1]->val[1] = %f\n", sl->snap[1]->val[1]);

//	for (n=0; n < SNAPLEN; n++)
//		sprintf (label, "%d", n+1);
//	sprintf (label, "\0");

	for (i=0; i < 1; i++) {						// find the mean of the noise - use only the first noise point
		for (j=0; j < noise->count; j++) {
			value = noise->snap[j]->raw[i] + 32767;		//make all values positive (range: -32767 to +32768)
			mean_noise += value / noise->count;		// calculate the mean for this point in the snapshot
		}
	}
	for (i=0; i < 1; i++) {					// calculate the standard deviation for the first point 
		tmp = 0;									//  in the noise snapshot
		sum = 0;
		for (j=0; j < noise->count; j++) {
			value = noise->snap[j]->raw[i] + 32767;		//make all values positive (range: -32767 to +32768)
			tmp = value - mean_noise, sum += tmp*tmp;
		}
		std_dev_noise = (double)(sqrt ((double)(sum / (noise->count-1))));
		binwidth_noise = 2*SD*std_dev_noise/100;
	}
	for (i=0; i < 1; i++) {
		lo_end_noise = mean_noise - (SD*std_dev_noise);
		//		printf ("noise: lo_end = %f; hi_end = %f\n", lo_end_noise, hi_end_noise);
//		PAUSE;
		for (j=0; j < noise->count; j++) {
			value = noise->snap[j]->raw[i] + 32767;		//make all values positive (range: -32767 to +32768)
			if(j<5) {
			  //				printf ("value = %f\n", value);
			  //				PAUSE;
			}

			bin = (int)(((value - lo_end_noise)/binwidth_noise)+1);
			if (bin < 0)
			  point_noise[0]++;
			else if (bin > 99)
			  point_noise[99]++;
			else
			  point_noise[bin]++;

//			  if (lo_end_noise <= value && value <= hi_end_noise) {
//				  bin = (int)(((value - lo_end_noise)/binwidth_noise)+1);
//				  (bin < 100 && bin >= 0) || DIE;
//				  point_noise[bin]++;
//			  }
//			  if (value < lo_end_noise)					// plot outliers in last bins
//				  point_noise[0]++;
//			  if (value > hi_end_noise)
//				  point_noise[99]++;

		}
		if (0)
		  for (j=0; j < 100; j++)
			printf ("point_noise[%d] = %d\n", j, point_noise[j]);
	}
//	PAUSE;
	//		printf ("sl->count = %d; noise->count = %d; noise scale = %f\n", sl->count, noise->count, (double)(sl->count)/(double)(noise->count));
	//		PAUSE;
	for (n=0; n < 100; n++) {						//normalize the point_noise array
		point_noise[n] = (int)(point_noise[n]*((double)(sl->count/(double)noise->count)));
	}

	for (i=0; i < SNAPLEN; i++) {
		for (j=0; j < sl->count; j++) {
			value = sl->snap[j]->val[i] + 32767;		//make all values positive (range: -32767 to +32768)
			mean[i] += value / sl->count;					// calculate the mean for this point in the snapshot
		}
	}
	for (i=0; i < SNAPLEN; i++) {					// calculate the standard deviation for each point in the snapshot
		tmp = 0;
		sum = 0;
		for (j=0; j < sl->count; j++) {
			value = sl->snap[j]->val[i] + 32767;		//make all values positive (range: -32767 to +32768)
			tmp = value - mean[i], sum += tmp*tmp;
		}
		std_dev[i] = (double)(sqrt ((double)(sum / (sl->count-1))));
		binwidth[i] = 2*SD*std_dev[i]/100;
	}
	for (i=0; i < SNAPLEN; i++) {
		for (j=0; j < sl->count; j++) {
			value = sl->snap[j]->val[i] + 32767;		//make all values positive (range: -32767 to +32768)
			lo_end = mean[i] - (SD*std_dev[i]);
			
			bin = (int)(((value - lo_end)/binwidth[i])+1);
			if (bin < 0)
			  point[i][0]++;
			else if (bin > 99)
			  point[i][99]++;
			else
			  point[i][bin]++;

//			  if (lo_end <= value && value <= hi_end) {
//				  bin = (int)(((value - lo_end)/binwidth[i])+1);
//				  (bin < 100 && bin >= 0) || DIE;
//				  point[i][bin]++;
//			  }
//			  if (value < lo_end)		// plot outliers in last bins
//				  point[i][0]++;
//			  if (value > hi_end)
//				  point[i][99]++;
		}
	}
	//	printf ("find the max value in array point\n");
	for (i=0; i < SNAPLEN; i++) {			// find the maximum bin in all the data histogram arrays
		for (j=0; j < 100; j++) {
			if (point[i][j] > maxmax) {
				maxmax = point[i][j];
				//				printf ("maxmax = %d\n", maxmax);
			}
		}
	}
	//	printf ("compute scale\n");
	scale = (float)(100)/maxmax;
	//	printf ("maxmax = %d; scale = %20.15f\n", maxmax, scale);
	//	printf ("scale array point\n");
	for (i = 0; i < SNAPLEN; i++) {
		for (j=0; j < 100; j++) {
			point[i][j] *= scale;
		}
	}
	for (j=0; j < 100; j++) {
		point_noise[j] *= scale;			// scale the noise data, too
		//		printf ("scaled point_noise[%d] = %d\n", j, point_noise[j]);
	}
//	PAUSE;
	
	ii=0;
	//	printf ("draw the histograms\n");
	for (n=0; n < 8; n++) {			// plot histograms as 3 standard deviations around the mean w/ overlaid noise plot
		for (m=0; m < 8 && ii < SNAPLEN; m++, ii++) {
			x = x0 + m*150;
			y = y0 + n*150;
//			pdftext (f, x0+50, y0+30, label[ii]);
			fprintf (f, "%d %d m %d %d l %d %d l s\n", x, y-100, x, y, x+100, y);		// draw the axes
			for (j=0; j < 100; j++){		// draw the histogram plots
				fprintf (f, "%d %d m %d %d l s\n", x+j+1, y, x+j+1, y-point[ii][j]);
			}
			fprintf (f, "1000. 0. 0. setrgbcolor\n");		// draw the noise overlay as a connected line (red)
			for (j=1; j < 100; j++) {					
				fprintf (f, " %d %d m %d %d l s\n",	x+j-1, y-point_noise[j-1], x+j, y-point_noise[j]);
			}
			fprintf (f, "0. setgray\n");
		}
	}
	end_pdf_page(f);
	return;
}
