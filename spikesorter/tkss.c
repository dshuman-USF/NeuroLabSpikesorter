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

/* tkss.c */

#include "nde.h"
#include <tk.h>
#include <tcl.h>
#include <string.h>
#include <limits.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "dmx.h"
#include <X11/Xutil.h>

int malloc_debug;

int
parse_long (ClientData clientData,
            Tcl_Interp *interp,
            Tk_Window tkwin,
            const char *value,
            char *widgRec,
            int offset)
{
  long *wrp = (long *)(widgRec + offset);
  *wrp = strtol (value, 0, 0);
  return TCL_OK;
}

const char *
print_long (ClientData clientData,
            Tk_Window tkwin,
            char *widgRec,
            int offset,
            Tcl_FreeProc **freeProcPtr)
{
  long *wrp = (long *)(widgRec + offset);
  char *s;
  asprintf (&s, "%ld", *wrp);
  *freeProcPtr = (Tcl_FreeProc *)free;
  return s;
}

Tk_CustomOption longint_option = {
  parse_long,
  print_long,
  0
};

typedef struct {
  unsigned short ccn;
  unsigned short spikestart;
  int segcnt;
  unsigned short rect[SNAPLEN-1][2];
  unsigned short raw[SNAPLEN];
  unsigned short y;
} Clus;

typedef struct {
  unsigned short v[SNAPLEN];
  XPoint *point;
  int member_count;
  int unit0;
  unsigned short **y_ye;
  XSegment *segment;
  int nsegments;
  unsigned short (**rect)[2];
  XRectangle *rectangle;
} Waves;

struct MapItem {struct MapItem *next; int *map;};
typedef struct MapItem MapItem;

typedef struct {short a, b; float d; unsigned short data;} Plane;
typedef struct {int count; XArc *arc;} ArcList;
typedef struct {Plane *plane; ArcList *arclist; XArc *arc; int narc; int used; int a, b;} Splot;

typedef struct {double mean; float uzmin; float uzmax; unsigned short ccnt;} WDT_Hdr;
typedef struct {WDT_Hdr *hdr; Clus **clus;} WDT;
typedef struct {int plane_count; int used; float d;} SPL_Hdr;
typedef struct {int offset, count; short a, b;} SPL_Idx;
typedef struct {SPL_Hdr *hdr; Plane **plane; int *narc; SPL_Idx *idx; FILE *f;} SPL;

typedef struct {double loc; int id; int uidx;} Entry;
typedef struct {
  int count; Entry *spike; char *filename; char *chan_path; int doublet_start;
  MapItem *first_map; MapItem *last_map; int map_count;
} Pos;
typedef struct {int bincnt, maxbin; double binsize; int *bin;} Hist;

int *unit_list, uidx_count;
int *id_to_uidx, id_count;

static int *merge_map;
static int *panel;
static int **member;
static int member_alloc;
static int *member_count;
static double *shift;
static double yscale_mult = 1;
static int shiftpxl;

static GC *gc_list;

static int *colormap;
static char *color_list[] = 
  {
    "#ff0000", "#ffdb00", "#92ff00", "#24b6ff", "#db00ff", "#ff00ff", "#0000ff", "#00ffff",
    "#924900", "#00ff00", "#7f7f7f", "#7f7f00"
  };

/*
static char *color_list[] = {
  "#ffffff",
  "#ff0000",			// 1 
  "#ffdb00",			// 2 
  "#92ff00",			// 3 
  "#24b6ff",			// 4 
  "#db00ff",			// 5 
  "#ff00ff",			// 6 
  "#0000ff",			// 7 
  "#00ffff",			// 8 
  "#ffffff",
  "#924900",			// 10 
  "#ffffff",
  "#00ff00",			// 12 
  "#ffffff",
  "#ffffff",
  "#7f7f7f",			// 15 
  "#ffff00"			// 16 
};
*/
static int color_count = sizeof color_list / sizeof *color_list;
static double sampling_frequency = SF;

typedef struct {XPoint pline[SNAPLEN]; int color; short sample[SNAPLEN]; } MW;
typedef struct {
  int goodmin;
  int goodcount;
  int hdrsize;
  int samples_per_file;
  int samples_per_swipe_time;
  int bytes_per_swipe;
  int offset_in_swipe;
  int channels_per_board;
  int chidx;
  MW *mw;
  int count;
  MW *ref;
  int refcount;
  int *unit;
  int max;
  int start;
  int next;
  int last;
  FILE *f;
  Pos *pos;
} MWaves;

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

void nullproc (ClientData x)
{
}


typedef struct {
  Tk_Window tkwin;		/* Window that embodies the square.  NULL
				 * means window has been deleted but
				 * widget record hasn't been cleaned up yet. */
  Display *display;		/* X's token for the window's display. */
  Tcl_Interp *interp;		/* Interpreter associated with widget. */
  Tcl_Command widgetCmd;	/* Token for square's widget command. */
  int run;
  double time;
  Tcl_TimerToken ttoken;
  long data;
  int type;
  XColor *color;
  /*
     * Information used when displaying widget:
     */

  int borderWidth;		/* Width of 3-D border around whole widget. */
  Tk_3DBorder bgBorder;	/* Used for drawing background. */
  Tk_3DBorder fgBorder;	/* For drawing square. */
  int relief;			/* Indicates whether window as a whole is
				 * raised, sunken, or flat. */
  GC gc;			/* Graphics context for copying from
				 * off-screen pixmap onto screen. */
  int doubleBuffer;		/* Non-zero means double-buffer redisplay
				 * with pixmap;  zero means draw straight
				 * onto the display. */
  int updatePending;		/* Non-zero means a call to SquareDisplay
				 * has already been scheduled. */
  double yscale;
} Square;

double cov_rep_root[SNAPLEN][SNAPLEN];
static double cvec[SNAPLEN];
double mean, sd;


/*
 * Information used for argv parsing.
 */

static Tk_ConfigSpec configSpecs[] = {
    {TK_CONFIG_BORDER, "-background", "background", "Background",
	"#d9d9d9", Tk_Offset(Square, bgBorder), TK_CONFIG_COLOR_ONLY},
    {TK_CONFIG_BORDER, "-background", "background", "Background",
	"white", Tk_Offset(Square, bgBorder), TK_CONFIG_MONO_ONLY},
    {TK_CONFIG_SYNONYM, "-bd", "borderWidth", (char *) NULL,
	(char *) NULL, 0, 0},
    {TK_CONFIG_SYNONYM, "-bg", "background", (char *) NULL,
	(char *) NULL, 0, 0},
    {TK_CONFIG_PIXELS, "-borderwidth", "borderWidth", "BorderWidth",
	"5", Tk_Offset(Square, borderWidth), 0},
    {TK_CONFIG_INT, "-dbl", "doubleBuffer", "DoubleBuffer",
	"1", Tk_Offset(Square, doubleBuffer), 0},
    {TK_CONFIG_CUSTOM, "-data", "data", "Data", "0", Tk_Offset(Square, data), 0, &longint_option},
    {TK_CONFIG_INT, "-run", "run", "Run", "0", Tk_Offset(Square, run), 0},
    {TK_CONFIG_DOUBLE, "-time", "time", "Time", "0", Tk_Offset(Square, time), 0},
    {TK_CONFIG_INT, "-ttoken", "ttoken", "Ttoken", "0", Tk_Offset(Square, ttoken), 0},
    {TK_CONFIG_INT, "-type", "type", "Type", "0", Tk_Offset(Square, type), 0},
    {TK_CONFIG_COLOR, "-color", "color", "Color", "black", Tk_Offset(Square, color), 0},
    {TK_CONFIG_SYNONYM, "-fg", "foreground", (char *) NULL,
	(char *) NULL, 0, 0},
    {TK_CONFIG_BORDER, "-foreground", "foreground", "Foreground",
	"#b03060", Tk_Offset(Square, fgBorder), TK_CONFIG_COLOR_ONLY},
    {TK_CONFIG_BORDER, "-foreground", "foreground", "Foreground",
	"black", Tk_Offset(Square, fgBorder), TK_CONFIG_MONO_ONLY},
    {TK_CONFIG_RELIEF, "-relief", "relief", "Relief",
	"raised", Tk_Offset(Square, relief), 0},
    {TK_CONFIG_DOUBLE, "-yscale", "yscale", "Yscale", "1", Tk_Offset(Square, yscale), 0},
    {TK_CONFIG_END, (char *) NULL, (char *) NULL, (char *) NULL,
	(char *) NULL, 0, 0}
};

/*
 * Forward declarations for procedures defined later in this file:
 */

int SquareCmd _ANSI_ARGS_((ClientData clientData, Tcl_Interp *interp, int argc, const char **argv));
static void SquareCmdDeletedProc _ANSI_ARGS_((ClientData clientData));
static int SquareConfigure _ANSI_ARGS_((Tcl_Interp *interp, Square *squarePtr, int argc, const char **argv, int flags));
static void SquareDestroy _ANSI_ARGS_((char *memPtr));
static void SquareDisplay _ANSI_ARGS_((ClientData clientData));
static void SquareEventProc _ANSI_ARGS_((ClientData clientData, XEvent *eventPtr));
static int SquareWidgetCmd _ANSI_ARGS_((ClientData clientData, Tcl_Interp *, int argc, const char **argv));
static void WaveDisplay (ClientData clientData);
static void WaveConfigure (ClientData clientData);
static void SplotDisplay (ClientData clientData);
static void SplotConfigure (ClientData clientData);
static void MWaveDisplay (ClientData clientData);
static void MWaveConfigure (ClientData clientData);

void (*paint[])(ClientData) = {SquareDisplay, WaveDisplay, SplotDisplay, MWaveDisplay};
void (*resize[])(ClientData) = {nullproc, WaveConfigure, SplotConfigure, MWaveConfigure};


static void get_next_mwave (ClientData clientData);

static int *
new_map (Tcl_Interp *interp, Pos *pos)
{
  MapItem *mapitem;
  char *s;

  TCALLOC (mapitem, 1);
  TMALLOC (mapitem->map, id_count);
  if (pos->first_map == 0)
    pos->last_map = pos->first_map = mapitem;
  else {
    pos->last_map->next = mapitem;
    pos->last_map = mapitem;
  }
  if (asprintf (&s, "%d", pos->map_count) == -1) exit (1);
  Tcl_SetVar (interp, "latest_merge_backup", s, TCL_GLOBAL_ONLY);
  free (s);
  pos->map_count++;
  return mapitem->map;
}

static int
match_map (Tcl_Interp *interp, Pos *pos)
{
  MapItem *mapitem;
  int n = 0, match = -1;
  char *s;

  for (n = 0, mapitem = pos->first_map; mapitem && match == -1; mapitem = mapitem->next, n++)
    if (memcmp (merge_map, mapitem->map, id_count * sizeof *merge_map) == 0)
      match = n;

  if (asprintf (&s, "%d", match) == -1) exit (1);
  Tcl_SetVar (interp, "current_merge_version", match < 0 ? "" : s, TCL_GLOBAL_ONLY);
  free (s);
  return match;
}

static void
read_versions (Tcl_Interp *interp, Pos *pos)
{
  int version, n;
  char *mrgv, *mrg;
  FILE *f;
  int *map;
  mrg = change_filetype_trim (pos->filename, ".mrg", 0);

  map = new_map (interp, pos);
  for (n = 0; n < id_count; n++)
    map[n] = n;
  for (version = 1; ; version++) {
    if (asprintf (&mrgv, "%s.~%d~", mrg, version) == -1) exit (1);
    if ((f = fopen (mrgv, "rb")) == 0)
      break;
    free (mrgv);
    map = new_map (interp, pos);
    fread (map, sizeof *map, id_count, f) == id_count || DIE;
    fclose (f);
  }
}

static void
write_mrg (Tcl_Interp *interp, Pos *pos)
{
  char *mrg;
  FILE *f;

  mrg = change_filetype_trim (pos->filename, ".mrg", 0);
  (f = fopen (mrg, "wb")) || DIE;
  fwrite (merge_map, sizeof *merge_map, id_count, f) == id_count || DIE;
  fclose (f);
}

static int
write_shift (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  char *fname;
  FILE *f;

  Pos *pos = (Pos *) atol (argv[1]);
  pos != 0 || DIE;
  fname = change_filetype_trim (pos->filename, ".shift", 0);
  (f = fopen (fname, "wb")) || DIE;
  fwrite (shift, sizeof *shift, id_count, f) == id_count || DIE;
  fclose (f);
  return TCL_OK;
}

static void
new_mrg (Tcl_Interp *interp, Pos *pos)
{
  char *mrg, *cmd;
  int *map;

  if (match_map (interp, pos) == -1) {
    mrg = change_filetype_trim (pos->filename, ".mrg", 0);
    if (asprintf (&cmd, "cp --force --backup=t %s %s", mrg, mrg) == -1) exit (1);
    if (system (cmd));
    free (cmd);
    map = new_map (interp, pos);
    memcpy (map, merge_map, id_count * sizeof *map);
  }
}

static void
norm_snap (float *v)
{
  int i;
  float fmean = mean;

  for (i = 0; i < SNAPLEN; i++)
    v[i] = (v[i] - fmean) / sd;
}

static void
my_whiten_snap (float *v)
{
  int i, j;
  float tmp[SNAPLEN];
  double sum;
  float fmean = mean;

  for (i = 0; i < SNAPLEN; i++) {
    sum = 0;
    for (j = 0; j < SNAPLEN; j++)
      sum += (v[j]-fmean) * cov_rep_root[j][i];
    tmp[i] = (float) sum;
  }
  memcpy (v, tmp, SNAPLEN * sizeof *tmp);
}

double
distance (float *a, float *b)
{
  double sum = 0, tmp;
  int n;

  for (n = 0; n < SNAPLEN; n++)
    tmp = a[n] - b[n], sum += tmp*tmp;
  return sqrt(sum);
}

static int
read_whitening (char *filename)
{
  FILE *f;

  
  if ((f = fopen (change_filetype_trim (filename, ".wht", 0), "rb")) == 0) {
    printf ("error opening %s\n", change_filetype_trim (filename, ".wht", 0));
    perror ("");
    return 0;
  }
  fread (&mean, sizeof mean, 1, f) == 1 || DIE;
  fread (cov_rep_root, sizeof cov_rep_root, 1, f) == 1 || DIE;
  fread (cvec, sizeof cvec, 1, f) == 1 || DIE;
  sd = sqrt (cvec[0]);
  fclose (f);
  return 1;
}

static void
print_distances (char *filename, short *v)
{
  FILE *f;
  int ccnt, n;
  static float (*c)[SNAPLEN];
  static float (*d)[SNAPLEN];
  float w[SNAPLEN];
  
  if (!read_whitening (filename))
    return;
  (f = fopen (change_filetype_trim (filename, ".ctr", 0), "rb")) || DIE;
  fread (&ccnt, sizeof ccnt, 1, f) == 1 || DIE;
  (ccnt > 0 && ccnt < 10000) || DIE;
  TREALLOC (c, ccnt);
  TREALLOC (d, ccnt);
  fread (c, ccnt * sizeof c[0], 1, f) == 1 || DIE;
  fclose (f);
  memcpy (d, c, ccnt * sizeof c[0]);
  for (n = 0; n < SNAPLEN; n++)
    w[n] = v[n];
  my_whiten_snap (w);
  for (n  = 0; n < ccnt; n++)
    printf (" %9d", n + 1);
  printf ("\n");
  for (n  = 0; n < ccnt; n++) {
    my_whiten_snap (c[n]);
    printf (" %9.6f", distance (w, &c[n][0]));
  }
  printf ("\n");
  for (n = 0; n < SNAPLEN; n++)
    w[n] = v[n];
  norm_snap (w);
  for (n  = 0; n < ccnt; n++) {
    norm_snap (d[n]);
    printf (" %9.6f", distance (w, &d[n][0]));
  }
  printf ("\n");
}

static void
MWaveDisplay (ClientData clientData)
{
  Square *squarePtr = (Square *) clientData;
  Tk_Window tkwin = squarePtr->tkwin;
  Pixmap pm = None;
  Drawable d;
  MWaves *mwaves;

  squarePtr->updatePending = 0;
  if (!Tk_IsMapped(tkwin))
    return;

  if (squarePtr->doubleBuffer) {
    pm = Tk_GetPixmap(Tk_Display(tkwin), Tk_WindowId(tkwin), Tk_Width(tkwin), Tk_Height(tkwin), DefaultDepthOfScreen(Tk_Screen(tkwin)));
    d = pm;
  }
  else d = Tk_WindowId(tkwin);

  Tk_Fill3DRectangle(tkwin, d, squarePtr->bgBorder, 0, 0, Tk_Width(tkwin), Tk_Height(tkwin), squarePtr->borderWidth, squarePtr->relief);

  mwaves = (MWaves *)squarePtr->data;

  if (mwaves) {
    int n, i;

    if (gc_list == 0) {
      TMALLOC (gc_list, color_count + 1);
      for (n = 0; n < color_count; n++)
	gc_list[n] = Tk_GCForColor (Tk_GetColor (squarePtr->interp, tkwin, color_list[n]), d);
      gc_list[n] = Tk_GCForColor (Tk_GetColor (squarePtr->interp, tkwin, "black"), d);
    }
    XSync (Tk_Display(tkwin), 0);
    if (1)
      for (n = 0; n < mwaves->refcount; n++)
	if (mwaves->ref[n].color != -1)
	  XDrawLines (Tk_Display(tkwin), d, gc_list[colormap[mwaves->ref[n].color]],
                      mwaves->ref[n].pline, SNAPLEN, CoordModeOrigin);
    for (n = 0; n < mwaves->count; n++) {
      i = (mwaves->start + n) % mwaves->max;
      XDrawLines (Tk_Display(tkwin), d, gc_list[colormap[mwaves->mw[i].color]],
                  mwaves->mw[i].pline, SNAPLEN, CoordModeOrigin);
    }

    if (0)
    {
      double xinc, yscale, ymid;
      int x1, y1, x2, y2;
      double mean = 1700.1;
      double thr = 7991.0;

      yscale = (double)(Tk_Height(tkwin) - 2 * squarePtr->borderWidth) / 65536;
      xinc = (double)(Tk_Width(tkwin) - 2 * squarePtr->borderWidth) / SNAPLEN;
      ymid = (double)(Tk_Height(tkwin) - 2 * squarePtr->borderWidth) / 2;

      x1 = 0;
      x2 = (SNAPLEN - 1) * xinc;
      y1 = y2 = squarePtr->borderWidth + floor (ymid - (mean + thr) * yscale + .5);
      XDrawLine (Tk_Display(tkwin), d, gc_list[color_count], x1, y1, x2, y2);

      x1 = 0;
      x2 = (SNAPLEN - 1) * xinc;
      y1 = y2 = squarePtr->borderWidth + floor (ymid - mean * yscale + .5);
      XDrawLine (Tk_Display(tkwin), d, gc_list[color_count], x1, y1, x2, y2);

      x1 = 0;
      x2 = (SNAPLEN - 1) * xinc;
      y1 = y2 = squarePtr->borderWidth + floor (ymid - (mean - thr) * yscale + .5);
      XDrawLine (Tk_Display(tkwin), d, gc_list[color_count], x1, y1, x2, y2);

      x1 = x2 = squarePtr->borderWidth + floor (9 * xinc + .5);
      y1 = 32767;
      y2 = -32768;
      XDrawLine (Tk_Display(tkwin), d, gc_list[color_count], x1, y1, x2, y2);

      x1 = x2 = squarePtr->borderWidth + floor (10 * xinc + .5);
      y1 = 32767;
      y2 = -32768;
      XDrawLine (Tk_Display(tkwin), d, gc_list[color_count], x1, y1, x2, y2);

    }

    static GC center_gc;

    if (center_gc == 0) {
      static XColor *white_color, *black_color;
      XGCValues values;
      white_color = Tk_GetColor (squarePtr->interp, tkwin, "white");
      Tk_GCForColor (white_color, d);
      black_color = Tk_GetColor (squarePtr->interp, tkwin, "black");
      values.line_style = LineDoubleDash;
      values.foreground =  black_color->pixel;
      values.background =  white_color->pixel;
      values.dashes = 4;
      center_gc = XCreateGC (Tk_Display(tkwin), d, GCLineStyle|GCDashList|GCBackground|GCForeground, &values);
    }
    
    double xinc = (double)(Tk_Width(tkwin) - 2 * squarePtr->borderWidth) / SNAPLEN;
    int x1, y1, x2, y2;
    x1 = squarePtr->borderWidth + floor (PRESAMPLES * xinc + .5); y1 = squarePtr->borderWidth;
    x2 = x1;                                                      y2 = Tk_Height(tkwin) - 1 - squarePtr->borderWidth;
    XDrawLine (Tk_Display(tkwin), d, center_gc, x1, y1, x2, y2);
    
  }
  if (squarePtr->doubleBuffer) {
    XCopyArea(Tk_Display(tkwin), pm, Tk_WindowId(tkwin), squarePtr->gc, 0, 0, (unsigned) Tk_Width(tkwin), (unsigned) Tk_Height(tkwin), 0, 0);
    Tk_FreePixmap(Tk_Display(tkwin), pm);
  }
}

static void
MWaveConfigure (ClientData clientData)
{
  Square *squarePtr = (Square *) clientData;
  Tk_Window tkwin = squarePtr->tkwin;
  MWaves *mwaves;

  mwaves = (MWaves *)squarePtr->data;

  if (mwaves) {
    double xinc, yscale, ymid;
    int n, i, j, whiten;
    float v[SNAPLEN];
    const char *s;

    xinc = (double)(Tk_Width(tkwin) - 2 * squarePtr->borderWidth) / SNAPLEN;
    ymid = (double)(Tk_Height(tkwin) - 2 * squarePtr->borderWidth) / 2;

    s = Tcl_GetVar (squarePtr->interp, "show_whitened", TCL_GLOBAL_ONLY);
    whiten = (s && strcmp (s, "1") == 0 && read_whitening (mwaves->pos->filename));
    {
      double wscale = cov_rep_root[0][0] ? 65536 * cov_rep_root[0][0] : 70;
      yscale = (double)(Tk_Height(tkwin) - 2 * squarePtr->borderWidth) / (whiten ? wscale : 65536) * yscale_mult;
    }

    for (n = 0; n < mwaves->count; n++) {
      j = (mwaves->start + n) % mwaves->max;
      for (i = 0; i < SNAPLEN; i++)
	v[i] = mwaves->mw[j].sample[i];
      if (whiten)
	my_whiten_snap (v);
      for (i = 0; i < SNAPLEN; i++) {
	mwaves->mw[j].pline[i].x = squarePtr->borderWidth + floor (i * xinc + .5);
	mwaves->mw[j].pline[i].y = squarePtr->borderWidth + floor (ymid - v[i] * yscale + .5);
      }
    }

    for (n = 0; n < mwaves->refcount; n++) {
      static char *buf;
      const char *s;
      int id = mwaves->unit[n];
      (id > 0 && id < id_count) || DIE;
      s = "";
      if (panel[id] >= 0) {
	if (asprintf (&buf, "cb%d", panel[id]) == -1) exit (1);
	s = Tcl_GetVar (squarePtr->interp, buf, TCL_GLOBAL_ONLY);
        free (buf);
      }
      if (merge_map[id] == id && s && strcmp (s, "1") == 0)
	mwaves->ref[n].color = id;
      else
	mwaves->ref[n].color = -1;
      for (i = 0; i < SNAPLEN; i++)
	v[i] = mwaves->ref[n].sample[i];
      if (whiten)
	my_whiten_snap (v);
      for (i = 0; i < SNAPLEN; i++) {
	mwaves->ref[n].pline[i].x = squarePtr->borderWidth + floor (i * xinc + .5);
	mwaves->ref[n].pline[i].y = squarePtr->borderWidth + floor (ymid - v[i] * yscale + .5);
      }
    }

    if (squarePtr->run) {
      if (!squarePtr->ttoken)
	get_next_mwave (clientData);
    }
    else {
      Tcl_DeleteTimerHandler (squarePtr->ttoken);
      squarePtr->ttoken = 0;
    }
  }
}

#define BUFSZ 20
static void
get_next_mwave (ClientData clientData)
{
  Square *squarePtr = (Square *) clientData;
  Tk_Window tkwin = squarePtr->tkwin;
  double xinc, yscale, ymid;
  int i, j, n = 0, search_count;
  MWaves *mwaves;
  static char buf[BUFSZ];
  int sample, swipe, samples_in_swipe, bytes_per_sample, sample_in_swipe_time;
  long long offset;
  static short *rbuf_tmp;
  short *rbuf;
  double binstart, binend;

  if (!squarePtr->run) {
    Tcl_DeleteTimerHandler (squarePtr->ttoken);
    squarePtr->ttoken = 0;
    return;
  }
  binstart = atoi (Tcl_GetVar (squarePtr->interp, "binstart", TCL_GLOBAL_ONLY));
  binend = atoi (Tcl_GetVar (squarePtr->interp, "binend", TCL_GLOBAL_ONLY));

  mwaves = (MWaves *)squarePtr->data;
  if (squarePtr->time == 0) {
    int mcnt, old_count, new_count, goal, fwd, new_last, inc;
    mcnt = atoi (Tcl_GetVar (squarePtr->interp, "movie_count", TCL_GLOBAL_ONLY));
    if (mcnt < 1)
      mcnt = 1;
    if (mcnt > mwaves->max)
      mcnt = mwaves->max;

    old_count = mwaves->count;
    new_count = mwaves->count < mcnt ? mwaves->count + 1 : mcnt;
    fwd = squarePtr->run > 0;
    
    if (mwaves->last == -1) {
      goal = fwd ? 0 : -2;
      n = mwaves->start;
      new_last = 0;
    }
    else if (fwd) {
      goal = old_count;
      mwaves->start = (mwaves->start + old_count - new_count + 1) % mwaves->max;
      n = (mwaves->start + new_count - 1)  % mwaves->max;
      new_last = new_count - 1;
    }
    else {
      goal =  MAX (old_count - new_count - 1, -1);
      mwaves->start = ((mwaves->start + goal) + mwaves->max) % mwaves->max;
      n = mwaves->start;
      new_last = 0;
    }
    inc = goal > mwaves->last ? 1 : -1;
    search_count = 0;
    for (i = mwaves->last; i != goal;) {
      const char *s;
      double mstime;
      int pnl;
      mwaves->next = mwaves->goodmin + (mwaves->next - mwaves->goodmin + inc + mwaves->goodcount) % mwaves->goodcount;
      double su = shift[mwaves->pos->spike[mwaves->next].id];
      mstime = (mwaves->pos->spike[mwaves->next].loc + su) / (SF / 1000);
      pnl = panel[mwaves->pos->spike[mwaves->next].id];
      s = "";
      if (pnl >= 0) {
	sprintf (buf, "cb%d", pnl);
	s = Tcl_GetVar (squarePtr->interp, buf, TCL_GLOBAL_ONLY);
      }
      if (strcmp (s, "1") == 0 && mstime >= binstart && mstime <= binend)
	i += inc;
      if (search_count++ > mwaves->goodcount) {
	squarePtr->run = 0;
	Tcl_DeleteTimerHandler (squarePtr->ttoken);
	squarePtr->ttoken = 0;
	return;
      }
    }
    mwaves->last = new_last;
    mwaves->count = new_count;
  }

  if (squarePtr->time != 0)
    sample = (int)floor (squarePtr->time * SF + .5);
  else {
    double su = shift[mwaves->pos->spike[mwaves->next].id];
    sample = floor ((mwaves->pos->spike[mwaves->next].loc + su) + .5) - PRESAMPLES;
  }
  if (sample >= mwaves->samples_per_file - SNAPLEN)
    sample = mwaves->samples_per_file - SNAPLEN;
  if (sample < 0)
    sample = 0;

  if (squarePtr->time != 0) {
    double d, min;
    int n_at_min = 0;
    min =  mwaves->samples_per_file;
    for (i = 0; i < mwaves->goodcount; i++) {
      n = mwaves->goodmin + i;
      double su = shift[mwaves->pos->spike[n].id];
      d = floor ((mwaves->pos->spike[n].loc + su) + .5) - PRESAMPLES - sample;
      if (fabs (d) < min)
	min = fabs (d),
	  n_at_min = n;
      if (d > 0)
	break;
    }
    mwaves->next = n_at_min;
    mwaves->count = 1;
    mwaves->last = -1;
    n = mwaves->start;
  }

  swipe = sample / mwaves->samples_per_swipe_time;
  sample_in_swipe_time = sample % mwaves->samples_per_swipe_time;
  samples_in_swipe = (sample_in_swipe_time + SNAPLEN > mwaves->samples_per_swipe_time
		      ? mwaves->samples_per_swipe_time - sample_in_swipe_time
		      : SNAPLEN);
  bytes_per_sample = 2;
  offset = (mwaves->hdrsize
	    + (long long) swipe * mwaves->bytes_per_swipe
	    + mwaves->offset_in_swipe
	    + sample_in_swipe_time * mwaves->channels_per_board * bytes_per_sample);
  fseeko64 (mwaves->f, offset, SEEK_SET);
  if (mwaves->channels_per_board > 1)
    TREALLOC (rbuf_tmp, SNAPLEN * mwaves->channels_per_board);
  rbuf = (mwaves->channels_per_board == 1) ? mwaves->mw[n].sample : rbuf_tmp;
  rbuf || DIE;
  samples_in_swipe <= SNAPLEN || DIE;
  fread64 (rbuf, mwaves->channels_per_board * bytes_per_sample, samples_in_swipe, mwaves->f) == samples_in_swipe || DIE;
  if (samples_in_swipe < SNAPLEN) {
    offset = mwaves->hdrsize + (long long)++swipe * mwaves->bytes_per_swipe + mwaves->offset_in_swipe;
    fseeko64 (mwaves->f, offset, SEEK_SET);
    fread64 (rbuf + samples_in_swipe * mwaves->channels_per_board,
	   mwaves->channels_per_board * bytes_per_sample,
	   SNAPLEN - samples_in_swipe, mwaves->f) == SNAPLEN - samples_in_swipe || DIE;
  }
  if (mwaves->channels_per_board > 1)
    for (i = 0, j = mwaves->chidx; i < SNAPLEN; i++, j += mwaves->channels_per_board)
      mwaves->mw[n].sample[i] = rbuf_tmp[j];

  if (0)
    if (fabs (squarePtr->run) > 1)
      print_distances (mwaves->pos->filename, mwaves->mw[n].sample);

  xinc = (double)(Tk_Width(tkwin) - 2 * squarePtr->borderWidth) / SNAPLEN;
  ymid = (double)(Tk_Height(tkwin) - 2 * squarePtr->borderWidth) / 2;

  {
    float v[SNAPLEN];
    const char *s;
    for (i = 0; i < SNAPLEN; i++)
      v[i] = mwaves->mw[n].sample[i];
    s = Tcl_GetVar (squarePtr->interp, "show_whitened", TCL_GLOBAL_ONLY);
    if (s && strcmp (s, "1") == 0 && read_whitening (mwaves->pos->filename)) {
      double wscale = cov_rep_root[0][0] ? 65536 * cov_rep_root[0][0] : 70;
      my_whiten_snap (v);
      yscale = (double)(Tk_Height(tkwin) - 2 * squarePtr->borderWidth) / wscale * yscale_mult;
    }
    else yscale = (double)(Tk_Height(tkwin) - 2 * squarePtr->borderWidth) / 65536 * yscale_mult;
    for (i = 0; i < SNAPLEN; i++) {
      mwaves->mw[n].pline[i].x = squarePtr->borderWidth + floor (i * xinc + .5);
      mwaves->mw[n].pline[i].y = squarePtr->borderWidth + floor (ymid - v[i] * yscale + .5);
    }
  }
  {
    int mm = merge_map[unit_list[mwaves->pos->spike[mwaves->next].uidx]];
    mm == -1 || mm > 0 || DIE;
    if (mm == -1)
      mm = 1;
    mwaves->mw[n].color = mm;
  }
  paint[squarePtr->type] (squarePtr);
  while (Tcl_DoOneEvent (TCL_ALL_EVENTS|TCL_DONT_WAIT))
    ;

  if (squarePtr->time == 0) {
    snprintf (buf, BUFSZ, "%d", mwaves->next + 1);
    Tcl_SetVar (squarePtr->interp, "spikenum", buf, TCL_GLOBAL_ONLY);
  }
  else {
    mwaves->mw[n].color = id_count;
    snprintf (buf, BUFSZ, "(%d)", mwaves->next + 1);
    Tcl_SetVar (squarePtr->interp, "spikenum", buf, TCL_GLOBAL_ONLY);
    squarePtr->time = 0;
  }

  snprintf (buf, BUFSZ, "%.6f", (double)sample/SF);
  Tcl_SetVar (squarePtr->interp, "spiketime", buf, TCL_GLOBAL_ONLY);

  snprintf (buf, BUFSZ, "%d",  mwaves->mw[n].sample[10]);
  Tcl_SetVar (squarePtr->interp, "spikeval", buf, TCL_GLOBAL_ONLY);

  {
    static int last_sample, last_run;
    if (squarePtr->run != last_run)
      last_run = squarePtr->run;
    else if ((squarePtr->run == 1 && sample < last_sample) || (squarePtr->run == -1 && sample > last_sample)) {
      squarePtr->run = 0;
      Tcl_SetVar (squarePtr->interp, "stopgo", "GO", TCL_GLOBAL_ONLY);
    }
    last_sample = sample;
  }

  if (abs (squarePtr->run) == 1)
    squarePtr->ttoken = Tcl_CreateTimerHandler(15, get_next_mwave, clientData);
  else {
    squarePtr->run = 0;
    Tcl_DeleteTimerHandler (squarePtr->ttoken);
    squarePtr->ttoken = 0;
  } 
  
}

static long
filesize (FILE *f)
{
  long bytes;

  fseek (f, 0, SEEK_END);
  bytes = ftell (f);
  rewind (f);
  return bytes;
}

static long long
filesize64 (FILE *f)
{
  long bytes;

  fseeko64 (f, 0, SEEK_END);
  bytes = ftello64 (f);
  fseeko64 (f, 0, SEEK_SET);
  return bytes;
}

static void
goodlimits (MWaves *mwaves, int samples_per_file)
{
  int goodmax;
  mwaves->goodmin = 0;
  double su = shift[mwaves->pos->spike[mwaves->goodmin].id];
  while (mwaves->goodmin < mwaves->pos->count && floor ((mwaves->pos->spike[mwaves->goodmin].loc + su) + .5) - PRESAMPLES < 0) {
    mwaves->goodmin++;
    su = shift[mwaves->pos->spike[mwaves->goodmin].id];
  }
  mwaves->goodmin < mwaves->pos->count || DIE;
  goodmax = mwaves->pos->count - 1;
  su = shift[mwaves->pos->spike[goodmax].id];
  while (goodmax >= 0 && floor ((mwaves->pos->spike[goodmax].loc + su) + .5) - PRESAMPLES + SNAPLEN > samples_per_file) {
     goodmax--;
     su = shift[mwaves->pos->spike[goodmax].id];
  }
  goodmax >= mwaves->goodmin || DIE;
  mwaves->goodcount = goodmax - mwaves->goodmin + 1;
}

static int
wave_params (char *file_name, MWaves *mwaves)
{
  char *fname, *dot, *chan_txt;
  FILE *f;
  static GlobalHeader *ghdr;
  static BoardHeader *bhdr;
  static int *enabled_channels, *swipe_size;
  int n, i, j, sample_rate, chan, chidx = -1;
  int board, channel, fpos;
  int swipe_size_total, bytes_per_sample;
  int skip, samples_per_swipe_time, bytes_per_swipe_time;

  bytes_per_sample = 2;
  if ((f = fopen64 (file_name, "rb"))) {
    long long fsize = filesize64 (f);
    mwaves->f = f;
    goodlimits (mwaves, fsize / bytes_per_sample);
    mwaves->samples_per_file = fsize / bytes_per_sample;
    mwaves->hdrsize = 0;
    mwaves->samples_per_swipe_time = fsize / bytes_per_sample;
    mwaves->bytes_per_swipe = fsize;
    mwaves->offset_in_swipe = 0;
    mwaves->channels_per_board = 1;
    mwaves->chidx = 0;
    return 1;
  }

  (fname = strdup (file_name)) || DIE;
  if (!(dot = strrchr (fname, '.')))
    dot = fname + strlen (fname);
  if (dot - fname < 3)
    goto fail;
  chan_txt = dot - 2;
  if (!isdigit ((int)chan_txt[0]) || !isdigit ((int)chan_txt[1]))
    goto fail;
  if ((chan = atoi (chan_txt)) == 0)
    goto fail;
  chan--;
  memmove (chan_txt, dot, strlen (dot) + 1);
  printf ("fname: %s\n", fname);
  if (!(f = fopen64 (fname, "rb")))
    goto fail;
  free (fname);
  mwaves->f = f;

  TREALLOC (ghdr, 1);
  if (fread64 (ghdr, sizeof *ghdr, 1, f) != 1
      || strncmp (ghdr->fid, "DATAMAX ", 8) != 0
      || chan >= ghdr->channel_count) {
    fclose64 (f);
    return 0;
  }
  board = chan / 8;
  (chan < ghdr->channel_count && chan >= 0) || DIE;
  channel = chan % 8;
  TREALLOC (bhdr, ghdr->board_count);
  TREALLOC (enabled_channels, ghdr->board_count);
  fread64 (bhdr, sizeof *bhdr, ghdr->board_count, f) == ghdr->board_count || DIE;
  for (sample_rate = n = 0; n < ghdr->board_count; n++) {
    enabled_channels[n] = 0;
    for (j = i = 0; i < 8; i++) {
      if (bhdr[n].chan[i].enable) {
	enabled_channels[n]++;
	if (sample_rate == 0)
	  sample_rate = bhdr[n].chan[i].sample_rate;
	else {
	  bhdr[n].chan[i].sample_rate == sample_rate || DIE;
	}
	bhdr[n].chan[i].board == n || DIE;
	bhdr[n].chan[i].channel == i || DIE;
	if (n == board && i == channel)
	  chidx = j;
	j++;
      }
    }
    enabled_channels[n] <= bhdr[n].channel_count || DIE;
  }
  chidx >= 0 || DIE;
  fpos = ftello64 (f);
  fpos == sizeof (*ghdr) + sizeof (*bhdr) * ghdr->board_count || DIE;
  fpos = (fpos / (1 << 16) + 1) * (1 << 16);

  TREALLOC (swipe_size, ghdr->board_count);
  samples_per_swipe_time = sample_rate / 10;
  bytes_per_swipe_time = samples_per_swipe_time * sizeof (short);
  for (swipe_size_total = n = 0; n < ghdr->board_count; n++)
    swipe_size_total += (swipe_size[n] = enabled_channels[n] * bytes_per_swipe_time);
  for (skip = n = 0; n < board; n++)
    skip += swipe_size[n];
  goodlimits (mwaves, ghdr->swipe_count * samples_per_swipe_time);
  mwaves->samples_per_file = ghdr->swipe_count * samples_per_swipe_time;
  mwaves->hdrsize = fpos;
  mwaves->samples_per_swipe_time = samples_per_swipe_time;
  mwaves->bytes_per_swipe = swipe_size_total;
  mwaves->offset_in_swipe = skip;
  mwaves->channels_per_board = enabled_channels[board];
  mwaves->chidx = chidx;

  return 1;

 fail:
  free (fname);
  return 0;
}


static int
get_mwaves(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  WDT *wdt;
  Pos *pos;
  MWaves *mwaves;
  int max, uidx;
  static char *buf;
  
  if (argc < 4)
    return TCL_OK;
  pos = (Pos *) atol (argv[1]);
  max = atoi (argv[2]);
  if ((mwaves = (MWaves *) atol (argv[3])) == 0)
    TCALLOC (mwaves, 1);
  wdt = (WDT *) atol (argv[4]);
  TREALLOC (mwaves->mw, max);
  TREALLOC (mwaves->ref, wdt->hdr->ccnt);
  TREALLOC (mwaves->unit, wdt->hdr->ccnt);
  for (uidx = 0; uidx < wdt->hdr->ccnt; uidx++) {
    int i;
    memcpy (mwaves->ref[uidx].sample, wdt->clus[uidx]->raw, SNAPLEN * sizeof mwaves->ref[uidx].sample[0]);
    mwaves->unit[uidx] = wdt->clus[uidx]->ccn + 1;
    for (i = 0; i < SNAPLEN; i++)
      mwaves->ref[uidx].sample[i] =  32767 - mwaves->ref[uidx].sample[i] / 1200.0 * 65535;
  }
  mwaves->refcount = wdt->hdr->ccnt;
  mwaves->max = max;
  mwaves->count = mwaves->start = 0;
  mwaves->last = -1;
  mwaves->next = -1;
  if (mwaves->f)
    fclose64 (mwaves->f);
  mwaves->pos = pos;
  wave_params (pos->chan_path, mwaves) || DIE;
  mwaves->next = mwaves->goodmin - 1;

  if (asprintf (&buf, "%ld", (long)mwaves) == -1) exit (1);
  Tcl_SetResult (interp, buf, (Tcl_FreeProc *)free);

  return TCL_OK;
}


/*
 *--------------------------------------------------------------
 *
 * SquareCmd --
 *
 *	This procedure is invoked to process the "square" Tcl
 *	command.  It creates a new "square" widget.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	A new widget is created and configured.
 *
 *--------------------------------------------------------------
 */

int
SquareCmd(clientData, interp, argc, argv)
    ClientData clientData;	/* Main window associated with
				 * interpreter. */
    Tcl_Interp *interp;		/* Current interpreter. */
    int argc;			/* Number of arguments. */
    const char **argv;		/* Argument strings. */
{
    Tk_Window mainwin = (Tk_Window) clientData;
    Square *squarePtr;
    Tk_Window tkwin;

    if (argc < 2) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
		argv[0], " pathName ?options?\"", (char *) NULL);
	return TCL_ERROR;
    }

    tkwin = Tk_CreateWindowFromPath(interp, mainwin, argv[1], (char *) NULL);

    if (tkwin == NULL) {
	return TCL_ERROR;
    }

    /*
     * Allocate and initialize the widget record.
     */

    squarePtr = (Square *) ckalloc(sizeof(Square));
    squarePtr->tkwin = tkwin;
    squarePtr->display = Tk_Display(tkwin);
    squarePtr->interp = interp;
    squarePtr->widgetCmd = Tcl_CreateCommand(interp,
	    Tk_PathName(squarePtr->tkwin), SquareWidgetCmd,
	    (ClientData) squarePtr, SquareCmdDeletedProc);
    squarePtr->borderWidth = 0;
    squarePtr->bgBorder = NULL;
    squarePtr->fgBorder = NULL;
    squarePtr->relief = TK_RELIEF_FLAT;
    squarePtr->gc = None;
    squarePtr->doubleBuffer = 1;
    squarePtr->updatePending = 0;
    squarePtr->color = 0;
    squarePtr->ttoken = 0;
    squarePtr->run = 0;
    squarePtr->time = 0;
    squarePtr->data = 0;
    squarePtr->type = 0;
    squarePtr->yscale = 0;

    Tk_CreateEventHandler(squarePtr->tkwin, ExposureMask|StructureNotifyMask,
	    SquareEventProc, (ClientData) squarePtr);
    if (SquareConfigure(interp, squarePtr, argc-2, argv+2, 0) != TCL_OK) {
	Tk_DestroyWindow(squarePtr->tkwin);
	return TCL_ERROR;
    }

    Tcl_SetResult (interp, Tk_PathName(squarePtr->tkwin), TCL_VOLATILE);
    return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * SquareWidgetCmd --
 *
 *	This procedure is invoked to process the Tcl command
 *	that corresponds to a widget managed by this module.
 *	See the user documentation for details on what it does.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	See the user documentation.
 *
 *--------------------------------------------------------------
 */

static int
SquareWidgetCmd(clientData, interp, argc, argv)
    ClientData clientData;		/* Information about square widget. */
    Tcl_Interp *interp;			/* Current interpreter. */
    int argc;				/* Number of arguments. */
    const char **argv;			/* Argument strings. */
{
    Square *squarePtr = (Square *) clientData;
    int result = TCL_OK;
    size_t length;
    char c;

    if (argc < 2) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
		argv[0], " option ?arg arg ...?\"", (char *) NULL);
	return TCL_ERROR;
    }
    Tcl_Preserve((ClientData) squarePtr);
    c = argv[1][0];
    length = strlen(argv[1]);
    if ((c == 'c') && (strncmp(argv[1], "cget", length) == 0)
	    && (length >= 2)) {
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
		    argv[0], " cget option\"",
		    (char *) NULL);
	    goto error;
	}
	result = Tk_ConfigureValue(interp, squarePtr->tkwin, configSpecs,
		(char *) squarePtr, argv[2], 0);
    } else if ((c == 'c') && (strncmp(argv[1], "configure", length) == 0)
	    && (length >= 2)) {
	if (argc == 2) {
	    result = Tk_ConfigureInfo(interp, squarePtr->tkwin, configSpecs,
		    (char *) squarePtr, (char *) NULL, 0);
	} else if (argc == 3) {
	    result = Tk_ConfigureInfo(interp, squarePtr->tkwin, configSpecs,
		    (char *) squarePtr, argv[2], 0);
	} else {
	    result = SquareConfigure(interp, squarePtr, argc-2, argv+2,
		    TK_CONFIG_ARGV_ONLY);
	}
    } else {
	Tcl_AppendResult(interp, "bad option \"", argv[1],
		"\": must be cget or configure",
		(char *) NULL);
	goto error;
    }
    if (!squarePtr->updatePending) {
	Tcl_DoWhenIdle(paint[squarePtr->type], (ClientData) squarePtr);
	squarePtr->updatePending = 1;
    }
    Tcl_Release((ClientData) squarePtr);
    return result;

    error:
    Tcl_Release((ClientData) squarePtr);
    return TCL_ERROR;
}

/*
 *----------------------------------------------------------------------
 *
 * SquareConfigure --
 *
 *	This procedure is called to process an argv/argc list in
 *	conjunction with the Tk option database to configure (or
 *	reconfigure) a square widget.
 *
 * Results:
 *	The return value is a standard Tcl result.  If TCL_ERROR is
 *	returned, then interp->result contains an error message.
 *
 * Side effects:
 *	Configuration information, such as colors, border width,
 *	etc. get set for squarePtr;  old resources get freed,
 *	if there were any.
 *
 *----------------------------------------------------------------------
 */

static int
SquareConfigure(interp, squarePtr, argc, argv, flags)
     Tcl_Interp *interp;			/* Used for error reporting. */
     Square *squarePtr;			/* Information about widget. */
     int argc;				/* Number of valid entries in argv. */
     const char **argv;			/* Arguments. */
     int flags;				/* Flags to pass to
					 * Tk_ConfigureWidget. */
{
  if (argc == 2 && strcmp (argv[0], "-data") == 0)
    squarePtr->data = atol (argv[1]);
  else
    if (Tk_ConfigureWidget(interp, squarePtr->tkwin, configSpecs,
                           argc, argv, (char *) squarePtr, flags) != TCL_OK) {
      return TCL_ERROR;
    }
  resize[squarePtr->type]((ClientData) squarePtr);

  /*
   * Set the background for the window and create a graphics context
   * for use during redisplay.
   */

  Tk_SetWindowBackground(squarePtr->tkwin,
                         Tk_3DBorderColor(squarePtr->bgBorder)->pixel);
  if (squarePtr->gc == None) {
    XGCValues gcValues;
    gcValues.function = GXcopy;
    gcValues.graphics_exposures = False;
    squarePtr->gc = Tk_GetGC(squarePtr->tkwin,
                             GCFunction|GCGraphicsExposures, &gcValues);
  }

  /*
   * Register the desired geometry for the window.  Then arrange for
   * the window to be redisplayed.
   */

  Tk_GeometryRequest(squarePtr->tkwin, 110, 10);
  Tk_SetInternalBorder(squarePtr->tkwin, squarePtr->borderWidth);
  if (!squarePtr->updatePending) {
    Tcl_DoWhenIdle(paint[squarePtr->type], (ClientData) squarePtr);
    squarePtr->updatePending = 1;
  }
  return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * SquareEventProc --
 *
 *	This procedure is invoked by the Tk dispatcher for various
 *	events on squares.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	When the window gets deleted, internal structures get
 *	cleaned up.  When it gets exposed, it is redisplayed.
 *
 *--------------------------------------------------------------
 */

static void
SquareEventProc(clientData, eventPtr)
    ClientData clientData;	/* Information about window. */
    XEvent *eventPtr;		/* Information about event. */
{
  Square *squarePtr = (Square *) clientData;

  if (eventPtr->type == Expose) {
    if (!squarePtr->updatePending) {
      Tcl_DoWhenIdle(paint[squarePtr->type], (ClientData) squarePtr);
      squarePtr->updatePending = 1;
    }
  }
  else if (eventPtr->type == ConfigureNotify) {
    resize[squarePtr->type]((ClientData) squarePtr);
    if (!squarePtr->updatePending) {
      Tcl_DoWhenIdle(paint[squarePtr->type], (ClientData) squarePtr);
      squarePtr->updatePending = 1;
    }
  }
  else if (eventPtr->type == DestroyNotify) {
    if (squarePtr->tkwin != NULL) {
      squarePtr->tkwin = NULL;
      Tcl_DeleteCommandFromToken(squarePtr->interp,
				 squarePtr->widgetCmd);
    }
    if (squarePtr->updatePending) {
      Tcl_CancelIdleCall(paint[squarePtr->type], (ClientData) squarePtr);
    }
    Tcl_EventuallyFree((ClientData) squarePtr, SquareDestroy);
  }
}

/*
 *----------------------------------------------------------------------
 *
 * SquareCmdDeletedProc --
 *
 *	This procedure is invoked when a widget command is deleted.  If
 *	the widget isn't already in the process of being destroyed,
 *	this command destroys it.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The widget is destroyed.
 *
 *----------------------------------------------------------------------
 */

static void
SquareCmdDeletedProc(clientData)
    ClientData clientData;	/* Pointer to widget record for widget. */
{
    Square *squarePtr = (Square *) clientData;
    Tk_Window tkwin = squarePtr->tkwin;
    /*
     * This procedure could be invoked either because the window was
     * destroyed and the command was then deleted (in which case tkwin
     * is NULL) or because the command was deleted, and then this procedure
     * destroys the widget.
     */

    if (tkwin != NULL) {
	squarePtr->tkwin = NULL;
	Tk_DestroyWindow(tkwin);
    }
}

/*
 *--------------------------------------------------------------
 *
 * SquareDisplay --
 *
 *	This procedure redraws the contents of a square window.
 *	It is invoked as a do-when-idle handler, so it only runs
 *	when there's nothing else for the application to do.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Information appears on the screen.
 *
 *--------------------------------------------------------------
 */

static void
SquareDisplay(clientData)
    ClientData clientData;	/* Information about window. */
{
    Square *squarePtr = (Square *) clientData;
    Tk_Window tkwin = squarePtr->tkwin;
    Pixmap pm = None;
    Drawable d;
    Hist *hist;

    squarePtr->updatePending = 0;
    if (!Tk_IsMapped(tkwin)) {
	return;
    }

    /*
     * Create a pixmap for double-buffering, if necessary.
     */

    if (squarePtr->doubleBuffer) {
	pm = Tk_GetPixmap(Tk_Display(tkwin), Tk_WindowId(tkwin),
		Tk_Width(tkwin), Tk_Height(tkwin),
		DefaultDepthOfScreen(Tk_Screen(tkwin)));
	d = pm;
    } else {
	d = Tk_WindowId(tkwin);
    }

    /*
     * Redraw the widget's background and border.
     */

    Tk_Fill3DRectangle(tkwin, d, squarePtr->bgBorder, 0, 0, Tk_Width(tkwin),
	    Tk_Height(tkwin), squarePtr->borderWidth, squarePtr->relief);

    /*
     * Display the square.
     */

    hist = (Hist *)squarePtr->data;
    if (hist) {
      int n, i;
      XPoint *point;
      double spikepixels, binpixels;
      static GC axis_gc;
      static XColor *axis_color;
      XGCValues values;

      if (hist->maxbin == 0)
	spikepixels = 0;
      else
	spikepixels = (Tk_Height(tkwin) - 2 * squarePtr->borderWidth) / (double)hist->maxbin;
      binpixels = (Tk_Width(tkwin) - 2 * squarePtr->borderWidth) / (double)hist->bincnt;
      //      printf ("spikepixels %f, maxbin %d\n", spikepixels, hist->maxbin);
      TMALLOC (point, hist->bincnt * 2 + 3);
      for (i = n = 0; n < hist->bincnt; n++, i += 2) {
	point[i].x = floor (n * binpixels + .5) + squarePtr->borderWidth;
	point[i+1].x = floor ((n + 1) * binpixels + .5) + squarePtr->borderWidth;
	point[i].y = point[i+1].y = Tk_Height(tkwin) - 1 - floor (hist->bin[n] * spikepixels + .5) - squarePtr->borderWidth;
	//	printf ("n %d, bin %d\n", n, hist->bin[n]);
      }
      values.foreground =  squarePtr->color->pixel;
      XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);
      XDrawLines (Tk_Display(tkwin), d, squarePtr->gc, point, hist->bincnt * 2, CoordModeOrigin);
      values.foreground =  BlackPixel (Tk_Display(tkwin), 0);
      XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);

      if (axis_gc == 0) {
	axis_color = Tk_GetColor (squarePtr->interp, tkwin, "red");
	axis_gc = Tk_GCForColor (axis_color, d);
      }
      point[0].x = squarePtr->borderWidth;                        point[0].y = squarePtr->borderWidth;
      point[1].x = squarePtr->borderWidth;                        point[1].y = Tk_Height(tkwin) - 1 - squarePtr->borderWidth;
      point[2].x = Tk_Width(tkwin) - 1 - squarePtr->borderWidth;  point[2].y = Tk_Height(tkwin) - 1 - squarePtr->borderWidth;
      XDrawLines (Tk_Display(tkwin), d, axis_gc, point, 3, CoordModeOrigin);
      free (point);
    }

    /*
     * If double-buffered, copy to the screen and release the pixmap.
     */

    if (squarePtr->doubleBuffer) {
	XCopyArea(Tk_Display(tkwin), pm, Tk_WindowId(tkwin), squarePtr->gc,
		0, 0, (unsigned) Tk_Width(tkwin), (unsigned) Tk_Height(tkwin),
		0, 0);
	Tk_FreePixmap(Tk_Display(tkwin), pm);
    }
}

static void
WaveDisplay (ClientData clientData)
{
  Square *squarePtr = (Square *) clientData;
  Tk_Window tkwin = squarePtr->tkwin;
  Pixmap pm = None;
  Drawable d;
  Waves *waves;

  squarePtr->updatePending = 0;
  if (!Tk_IsMapped(tkwin))
    return;

  if (squarePtr->doubleBuffer) {
    pm = Tk_GetPixmap(Tk_Display(tkwin), Tk_WindowId(tkwin), Tk_Width(tkwin), Tk_Height(tkwin), DefaultDepthOfScreen(Tk_Screen(tkwin)));
    d = pm;
  }
  else d = Tk_WindowId(tkwin);

  Tk_Fill3DRectangle(tkwin, d, squarePtr->bgBorder, 0, 0, Tk_Width(tkwin), Tk_Height(tkwin), squarePtr->borderWidth, squarePtr->relief);

  waves = (Waves *)squarePtr->data;

  if (waves) {
    XPoint *point;
    static GC axis_gc, center_gc;
    static XColor *axis_color, *white_color, *black_color;
    XGCValues values;


    if (axis_gc == 0) {
      axis_color = Tk_GetColor (squarePtr->interp, tkwin, "red");
      axis_gc = Tk_GCForColor (axis_color, d);
      white_color = Tk_GetColor (squarePtr->interp, tkwin, "white");
      Tk_GCForColor (white_color, d);
      black_color = Tk_GetColor (squarePtr->interp, tkwin, "black");
      values.line_style = LineDoubleDash;
      values.foreground =  black_color->pixel;
      values.background =  white_color->pixel;
      values.dashes = 4;
      center_gc = XCreateGC (Tk_Display(tkwin), d, GCLineStyle|GCDashList|GCBackground|GCForeground, &values);
    }

    values.foreground =  squarePtr->color->pixel;
    XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);
    if (0)
    {
      int n;
      XSegment *s = waves->segment;
      for (n = 0; n < waves->nsegments; n++)
	if (s[n].x1 == s[n].x2 && s[n].x1 == s[n].y1 && s[n].x1 == s[n].y2)
	  printf ("all same\n");
    }
    if (waves->segment) {
      int count = waves->nsegments;
      int start = 0;
      while (count > 500000) {
	XDrawSegments (Tk_Display(tkwin), d, squarePtr->gc, waves->segment + start, 500000);
	start += 500000;
	count -= 500000;
      }
      XDrawSegments (Tk_Display(tkwin), d, squarePtr->gc, waves->segment + start, count);
    }

    if (waves->rectangle)
      XFillRectangles (Tk_Display(tkwin), d, squarePtr->gc, waves->rectangle, waves->member_count * (SNAPLEN - 1));
    
    values.foreground =  black_color->pixel;
    XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);

    XDrawLines (Tk_Display(tkwin), d, center_gc, waves->point, SNAPLEN, CoordModeOrigin);

    TMALLOC (point, 3);
    point[0].x = squarePtr->borderWidth;                        point[0].y = squarePtr->borderWidth;
    point[1].x = squarePtr->borderWidth;                        point[1].y = Tk_Height(tkwin) - 1 - squarePtr->borderWidth;
    point[2].x = Tk_Width(tkwin) - 1 - squarePtr->borderWidth;  point[2].y = Tk_Height(tkwin) - 1 - squarePtr->borderWidth;

    XDrawLines (Tk_Display(tkwin), d, axis_gc, point, 3, CoordModeOrigin);

    double xinc = (double)(Tk_Width(tkwin) - 2 * squarePtr->borderWidth) / SNAPLEN;
    point[0].x = squarePtr->borderWidth + floor (20 * xinc + .5); point[0].y = squarePtr->borderWidth;
    point[1].x = point[0].x;                                      point[1].y = Tk_Height(tkwin) - 1 - squarePtr->borderWidth;
    XDrawLines (Tk_Display(tkwin), d, center_gc, point, 2, CoordModeOrigin);
    

    free (point);
  }

  if (squarePtr->doubleBuffer) {
    XCopyArea(Tk_Display(tkwin), pm, Tk_WindowId(tkwin), squarePtr->gc, 0, 0, (unsigned) Tk_Width(tkwin), (unsigned) Tk_Height(tkwin), 0, 0);
    Tk_FreePixmap(Tk_Display(tkwin), pm);
  }
}

static void
WaveConfigure (ClientData clientData)
{
  Square *squarePtr = (Square *) clientData;
  Tk_Window tkwin = squarePtr->tkwin;
  Waves *waves;

  waves = (Waves *)squarePtr->data;

  if (waves) {
    double xinc, yscale, ymid;
    int n, flagval, x, sidx, x1, x2, y1, y2, midx;
    unsigned short *buf;
    XSegment *seg;

    if (!waves->point)
      TMALLOC (waves->point, SNAPLEN);

    xinc = (double)(Tk_Width(tkwin) - 2 * squarePtr->borderWidth) / SNAPLEN;
    yscale = (double)(Tk_Height(tkwin) - 2 * squarePtr->borderWidth) / 1200 * yscale_mult;
    ymid = (double)(Tk_Height(tkwin) - 2 * squarePtr->borderWidth) / 2;

    double dshift = shiftpxl / xinc;
    for (midx = 0; midx < waves->member_count; midx++)
      shift[member[waves->unit0][midx]] += dshift;
    
    for (n = 0; n < SNAPLEN; n++) {
      waves->point[n].x = squarePtr->borderWidth + floor ((n - shift[waves->unit0]) * xinc + .5);
      waves->point[n].y = squarePtr->borderWidth + floor (((double)waves->v[n] - 600) * yscale + ymid + .5);
    }
    TREALLOC (waves->rectangle, waves->member_count * (SNAPLEN-1));
    sidx = 0;
    for (midx = 0; midx < waves->member_count; midx++) {
      double su = shift[member[waves->unit0][midx]];
      for (n = 0; n < SNAPLEN - 1; n++) {
	int i = (midx * (SNAPLEN-1)) + n;
	int r0 = waves->rect[midx][n][0];
	int r1 = waves->rect[midx][n][1];
	waves->rectangle[i].x = squarePtr->borderWidth + floor ((n - su) * xinc + .5);
	waves->rectangle[i].y = squarePtr->borderWidth + floor (((double)r0 - 600) * yscale + ymid + .5);
	waves->rectangle[i].width = floor (xinc + .5) + 1;
	if      (r1 <  r0) waves->rectangle[i].height = 0;
	else if (r1 == r0) waves->rectangle[i].height = 1;
	else               waves->rectangle[i].height = floor ((r1 - r0) * yscale + .5);
      }
      TREALLOC (waves->segment, waves->nsegments);
      buf = waves->y_ye[midx];
      seg = waves->segment;
      flagval = 0x8000;
      n = x = 0;
      while ((buf[n] & 0xC000) == flagval) {
	x1 = floor ((x++ - su) * xinc + .5) + squarePtr->borderWidth;
	x2 = floor ((x   - su) * xinc + .5) + squarePtr->borderWidth;
	while ((buf[n] & 0xC000) == flagval) {
	  flagval = 0x4000;
	  y1 = floor (((buf[n++] & 0x3FFF) - 600) * yscale + ymid + .5) + squarePtr->borderWidth;
	  while ((buf[n] & 0xC000) == 0) {
	    if (sidx >= waves->nsegments) {
	      printf ("line %d: %d segments\n", __LINE__, waves->nsegments);
	      exit (0);
	    }
	    seg[sidx].x1 = x1;
	    seg[sidx].y1 = y1;
	    seg[sidx].x2 = x2;
	    y2 = floor ((buf[n++] - 600) * yscale + ymid + .5) + squarePtr->borderWidth;
	    seg[sidx++].y2 = y2;
	  }
	}
	flagval = 0x8000;
      }
    }
  }
}

static void
SplotDisplay (ClientData clientData)
{
  Square *squarePtr = (Square *) clientData;
  Tk_Window tkwin = squarePtr->tkwin;
  Pixmap pm = None;
  Drawable d;
  Splot *splot;
  int n;
  static unsigned long color[sizeof color_list / sizeof *color_list], black, white;

  squarePtr->updatePending = 0;
  if (!Tk_IsMapped(tkwin))
    return;
  if (Tk_Height (tkwin) != Tk_Width(tkwin) / 4 * 3)
    return;
  if (Tk_Width (tkwin) != Tk_Width(Tk_Parent(tkwin)))
    return;

  if (squarePtr->doubleBuffer) {
    pm = Tk_GetPixmap(Tk_Display(tkwin), Tk_WindowId(tkwin), Tk_Width(tkwin), Tk_Height(tkwin), DefaultDepthOfScreen(Tk_Screen(tkwin)));
    d = pm;
  }
  else d = Tk_WindowId(tkwin);

  Tk_Fill3DRectangle(tkwin, d, squarePtr->bgBorder, 0, 0, Tk_Width(tkwin), Tk_Height(tkwin), squarePtr->borderWidth, squarePtr->relief);

  if (color[0] == 0) {
    for (n = 0; n < color_count; n++)
      color[n] = Tk_GetColor (squarePtr->interp, tkwin, color_list[n])->pixel;
    black = Tk_GetColor (squarePtr->interp, tkwin, "black")->pixel;
    white = Tk_GetColor (squarePtr->interp, tkwin, "white")->pixel;
  }

  splot = (Splot *)squarePtr->data;

  merge_map || DIE;
  if (splot) {
    int uidx;
    XGCValues values;

    if (!(splot->used == uidx_count + 1))
      printf ("splot->used: %d, uidx_count: %d\n", splot->used, uidx_count);
    if (0) {
      uidx = splot->used - 2;	/* unclassified */
      values.foreground = color[(unit_list[uidx] - 1) % color_count];
      XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);
      XFillArcs (Tk_Display(tkwin), d, squarePtr->gc, splot->arclist[uidx].arc,  splot->arclist[uidx].count);
    }

    if (0)
    for (uidx = 0; uidx < splot->used - 2; uidx++) {
      if (unit_list[uidx] == splot->a + 1 || unit_list[uidx] == splot->b + 1)
	continue;
      values.foreground =  color[(unit_list[uidx] - 1) % color_count];
      XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);
      XFillArcs (Tk_Display(tkwin), d, squarePtr->gc, splot->arclist[uidx].arc,  splot->arclist[uidx].count);
      values.foreground = black;
      XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);
      XFillArcs (Tk_Display(tkwin), d, squarePtr->gc, splot->arclist[uidx].arc + splot->arclist[uidx].count - 1, 1);
      if (0) {
	XArc *arc = splot->arclist[uidx].arc + splot->arclist[uidx].count - 1;
	printf ("%s line %d, uidx %d: %d %d\n", __FILE__, __LINE__, uidx, arc->x, arc->y);
      }
    }

    uidx = splot->used - 1;
    values.foreground = black;
    XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);
    XFillArcs (Tk_Display(tkwin), d, squarePtr->gc, splot->arclist[uidx].arc,  splot->arclist[uidx].count);

    values.foreground = white;
    XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);
    XFillArcs (Tk_Display(tkwin), d, squarePtr->gc, splot->arclist[uidx].arc + splot->arclist[uidx].count - 1, 1);
    if (0) {
      XArc *arc = splot->arclist[uidx].arc + splot->arclist[uidx].count - 1;
      printf ("%s line %d, uidx %d: %d %d\n", __FILE__, __LINE__, uidx, arc->x, arc->y);
    }

    for (uidx = 0; uidx < uidx_count; uidx++) {
      int ab;
      if      (merge_map[unit_list[uidx]] == splot->a + 1) ab = splot->a;
      else if (merge_map[unit_list[uidx]] == splot->b + 1) ab = splot->b;
      else continue;
      values.foreground =  color[colormap[ab + 1]];
      XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);
      XFillArcs (Tk_Display(tkwin), d, squarePtr->gc, splot->arclist[uidx].arc,  splot->arclist[uidx].count);
      values.foreground = black;
      XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);
      if (unit_list[uidx] == ab + 1)
	XFillArcs (Tk_Display(tkwin), d, squarePtr->gc, splot->arclist[uidx].arc + splot->arclist[uidx].count - 1, 1);
    }

    if (0)
      for (uidx = 0; uidx < splot->used - 2; uidx++) {
	if (unit_list[uidx] != 5)
	  continue;
	values.foreground =  color[(unit_list[uidx] - 1) % color_count];
	XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);
	XFillArcs (Tk_Display(tkwin), d, squarePtr->gc, splot->arclist[uidx].arc,  splot->arclist[uidx].count);
	values.foreground = black;
	XChangeGC (Tk_Display(tkwin), squarePtr->gc, GCForeground, &values);
	XFillArcs (Tk_Display(tkwin), d, squarePtr->gc, splot->arclist[uidx].arc + splot->arclist[uidx].count - 1, 1);
      }

  }

  if (squarePtr->doubleBuffer) {
    XCopyArea(Tk_Display(tkwin), pm, Tk_WindowId(tkwin), squarePtr->gc, 0, 0, (unsigned) Tk_Width(tkwin), (unsigned) Tk_Height(tkwin), 0, 0);
    Tk_FreePixmap(Tk_Display(tkwin), pm);
  }
}

static int
get_distance (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Splot *splot;
  static char *buf;

  if (argc < 2)
    return TCL_OK;
  if ((splot = (Splot *) atol (argv[1])) == 0)
    return TCL_OK;

  if (asprintf (&buf, "%.1f", splot->plane->d) == -1) exit (1);
  Tcl_SetResult (interp, buf, (Tcl_FreeProc *)free);
  
  return TCL_OK;
}


static void
SplotConfigure (ClientData clientData)
{
  Square *squarePtr = (Square *) clientData;
  Tk_Window tkwin = squarePtr->tkwin;
  Splot *splot;

  splot = (Splot *)squarePtr->data;

  if (splot) {
    double xscale, yscale;
    int bufidx, bufidx0, arcidx, uidx;
    unsigned short *buf;
    XArc *arc;

    xscale = (double)(Tk_Width(tkwin) - 2 * squarePtr->borderWidth) / 1600;
    yscale = (double)(Tk_Height(tkwin) - 2 * squarePtr->borderWidth) / 1200;
    buf = &splot->plane->data;
    arc = splot->arc;
    for (arcidx = bufidx = uidx = 0; uidx < splot->used; uidx++) {
      int ccn = buf[bufidx];
      splot->arclist[uidx].arc = &arc[arcidx];
      bufidx++;
      bufidx0 = bufidx;
      while (buf[bufidx] != 0x8000) {
	arc[arcidx].x = squarePtr->borderWidth + floor (buf[bufidx++] * xscale + .5);
	arc[arcidx].y = squarePtr->borderWidth + floor (buf[bufidx++] * yscale + .5);
	arcidx++;
      }
      if (0)
	printf ("%s line %d: ccn %d, x %d, y %d, x %d, y %d, mark %x\n", __FILE__, __LINE__, ccn,
		buf[bufidx - 4], buf[bufidx - 3],
		buf[bufidx - 2], buf[bufidx - 1],
		buf[bufidx]
		);
      splot->arclist[uidx].count = (bufidx - bufidx0) / 2;
      bufidx++;
    }
  }
}


/*
 *----------------------------------------------------------------------
 *
 * SquareDestroy --
 *
 *	This procedure is invoked by Tcl_EventuallyFree or Tcl_Release
 *	to clean up the internal structure of a square at a safe time
 *	(when no-one is using it anymore).
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Everything associated with the square is freed up.
 *
 *----------------------------------------------------------------------
 */

static void
SquareDestroy(memPtr)
    char *memPtr;		/* Info about square widget. */
{
    Square *squarePtr = (Square *) memPtr;
    if (squarePtr->ttoken) {
      Tcl_DeleteTimerHandler (squarePtr->ttoken);
      squarePtr->ttoken = 0;
      squarePtr->run = 0;
      return;
    }
    XSync (squarePtr->display, 0);

    Tk_FreeOptions(configSpecs, (char *) squarePtr, squarePtr->display, 0);
    if (squarePtr->gc != None) {
	Tk_FreeGC(squarePtr->display, squarePtr->gc);
    }
    ckfree((char *) squarePtr);
}

static int
PutIntCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  long val;

  if (argc > 1) {
    val = atoi (argv[1]);
    printf ("val: %ld\n", val);
    ((void (*)())val)();
  }
  else
    printf ("one object\n");
  return TCL_OK;
}

static int
get_units(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Pos *pos;
  int n, count, i, j, id;
  char *buf;
  FILE *f;

  if (argc < 2)
    return TCL_OK;
  pos = (Pos *) atol (argv[1]);
  for (count = n = 0; n < pos->count; n++) {
    for (i = 0; i < count; i++)
      if (unit_list[i] >= pos->spike[n].id)
	break;
    if (i < count && (unit_list[i] == pos->spike[n].id))
      continue;
    TREALLOC (unit_list, ++count);
    for (j = count - 1; j > i; j--)
      unit_list[j] = unit_list[j-1];
    unit_list[i] =  pos->spike[n].id;
  }
  uidx_count = count;
  if (unit_list == 0)
    return TCL_OK;
  unit_list[0] || DIE;
  id_count = unit_list[count - 1] + 1;
  TREALLOC (id_to_uidx, id_count);
  for (n = 0; n < count; n++)
    id_to_uidx[unit_list[n]] = n;
  for (count = n = 0; n < pos->count; n++)
    pos->spike[n].uidx = id_to_uidx[pos->spike[n].id];

  if (shift != NULL) TREALLOC (shift, id_count);
  else {
    TCALLOC (shift, id_count);
    FILE *f = fopen (change_filetype_trim (pos->filename, ".shift", 0), "rb");
    if (f != NULL) {
      fread (shift, sizeof *shift, id_count, f) == id_count || DIE;
      fclose (f);
    }
  }

  TREALLOC (merge_map, id_count);
  TREALLOC (panel, id_count);
  for (n = 0; n < member_alloc; n++)
    free (member[n]);
  TREALLOC (member, member_alloc = id_count);
  TMEMSET (member, 0, id_count);
  TREALLOC (member_count, id_count);
  if ((f = fopen (change_filetype_trim (pos->filename, ".mrg", 0), "rb")) == 0) {
    for (n = 0; n < id_count; n++)
      merge_map[n] = n;
  }
  else {
    fread (merge_map, sizeof *merge_map, id_count, f) == id_count || DIE;
    fclose (f);
  }

  TREALLOC (colormap, id_count + 1);
  for (int i = 0; i < id_count; ++i)
    colormap[i] = i % color_count;
  colormap[id_count] = color_count;

  if (pos->first_map == 0)
    read_versions (interp, pos);
  match_map (interp, pos);
  for (id = 0; id < id_count; id++) {
    member_count[id] = 0;
    for (n = 0; n < uidx_count; n++)
      if (merge_map[unit_list[n]] == id) {
	int midx = member_count[id]++;
	TREALLOC (member[id], member_count[id]);
	member[id][midx] = unit_list[n];
      }
  }
  for (n = 0; n < id_count; n++)
    panel[n] = -1;
  for (n = i = 0; n < uidx_count; n++) {
    id = unit_list[n];
    if (member_count[id]) {
      for (j = 0; j < member_count[id]; j++)
	panel[member[id][j]] = i;
      i++;
      if (asprintf (&buf, "%d", id) == -1) exit (1);
      Tcl_AppendElement (interp, buf);
      free (buf);
    }
  }
  return TCL_OK;
}

static int
write_edt(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Pos *pos;
  int mapsize, n, old, new, seen_count, i, maxmap;
  int *map;
  char *seen;
  FILE *f;

  if (argc < 4)
    return TCL_OK;
  (pos = (Pos *) atol (argv[1])) != 0 || DIE;
  (f = fopen (change_filetype_trim (argv[2], ".edt", 0), "w")) || DIE;
  fprintf (f, "   33   3333333\n");
  fprintf (f, "   33   3333333\n");
  mapsize =  atoi (argv[3]) + 1;
  TCALLOC (map, mapsize);
  TCALLOC (seen, mapsize);
  maxmap = 0;
  for (n = 4; n < argc; n += 2) {
    old = atoi (argv[n]);
    new = atoi (argv[n+1]);
    old < mapsize || DIE;
    map[old] = new;
    if (new > maxmap)
      maxmap = new;
  }
  typedef struct {int id; int t;} IDT;
  IDT *idt;
  TMALLOC (idt, pos->count);
  int ii = 0;
  for (seen_count = n = 0; n < pos->count; n++) {
    int m_id;
    (pos->spike[n].id > 0 && pos->spike[n].id < id_count) || DIE;
    (m_id = merge_map[pos->spike[n].id]) < mapsize || DIE;
    if (m_id > 0 && (new = map[m_id]) > 0) {
      if (!seen[m_id]) {
	seen[m_id] = 1;
	seen_count++;
      }
      double su = shift[pos->spike[n].id];
      idt[ii].id = new;
      idt[ii].t = (int)floor((pos->spike[n].loc + su)/(SF/10000)+.5);
      ++ii;
    }
  }
  int idtcnt = ii;
  int idtcmp (const void *va, const void *vb)
  {
    const IDT *a = (const IDT *) va;
    const IDT *b = (const IDT *) vb;

    if (a->id == b->id && a->t == b->t) return 0;
    if (a->t < b->t || (a->t == b->t && a->id < b->id)) return -1;
    return 1;
  }
  qsort (idt, idtcnt, sizeof *idt, idtcmp);
  for (n = 0; n < idtcnt; n++)
    fprintf (f, "%5d%10d\n", idt[n].id, idt[n].t);
  fclose (f);
  free (idt);
  (f = fopen (change_filetype_trim (argv[2], ".unt", 0), "w")) || DIE;
  fprintf (f, "%d\n", seen_count);
  for (i = 0; i <= maxmap; i++)
    for (n = 0; n < mapsize; n++)
      if (seen[n] && map[n] == i) {
	fprintf (f, "%d\n", i);
	break;
      }
  fclose (f);

  (f = fopen (change_filetype_trim (argv[2], ".status", 0), "a")) || DIE;
  putc ('E', f);
  fclose (f);

  free (map);
  free (seen);
  return TCL_OK;
}

static int
get_maxbin(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Hist *hist;
  char *buf;

  if (argc < 2)
    return TCL_OK;
  hist = (Hist *) atol (argv[1]);
  if (asprintf (&buf, "%d", hist->maxbin) == -1) exit (1);
  Tcl_AppendElement (interp, buf);
  free (buf);
  return TCL_OK;
}

static int
get_last_time (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Pos *pos;
  char *buf;
  double samples_per_ms = sampling_frequency / 1000;

  if (argc < 2)
    return TCL_OK;
  pos = (Pos *) atol (argv[1]);
  double su = shift ? shift[pos->spike[pos->count-1].id] : 0;
  if (asprintf (&buf, "%f", (pos->spike[pos->count-1].loc + su) / samples_per_ms) == -1) exit (1);
  Tcl_AppendElement (interp, buf);
  free (buf);
  
  return TCL_OK;
}

static int
wait_pid (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Tcl_Pid pid, retval;
  int options, stat;
  char *buf;

  if (argc < 2)
    return TCL_OK;
  pid = (Tcl_Pid) atol (argv[1]);
  options = 0;
  if (argc > 2 && strcmp (argv[2], "nohang") == 0)
    options = WNOHANG;

  retval = Tcl_WaitPid (pid, &stat, options);

  if (asprintf (&buf, "%ld", (long)retval) == -1) exit (1);
  Tcl_AppendElement (interp, buf);
  free (buf);
  
  if (asprintf (&buf, "%d", (int)stat) == -1) exit (1);
  Tcl_AppendElement (interp, buf);
  free (buf);
  
  return TCL_OK;
}

static int
get_color(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  int n;
  if (argc < 2)
    return TCL_OK;
  n = atoi (argv[1]);

  Tcl_SetResult (interp, color_list[colormap[n + 1]], TCL_STATIC);

  return TCL_OK;
}

static int
recolor (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  static gsl_rng *rng;
  if (rng == NULL) rng = gsl_rng_alloc (gsl_rng_ranlxs2);

  int count = 0;
  for (int i = 1; i < id_count; ++i)
    if (merge_map[i] == i)
      ++count;

  int setcnt = ceil (count / (1.0 * color_count));
  int *sets;
  TMALLOC (sets, setcnt * color_count);
  for (int iset = 0; iset < setcnt; ++iset)
    for (int i = 0; i < color_count; ++i)
      sets[iset * color_count + i] = i;

  for (int iset = 0; iset < setcnt; ++iset)
    gsl_ran_shuffle (rng, sets + iset * color_count, color_count, sizeof *sets);

  int setsi = 0;
  for (int i = 1; i < id_count; ++i)
    if (merge_map[i] == i)
      colormap[i] = sets[setsi++];
  free (sets);
  
  return TCL_OK;
}

static int
get_waves(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  WDT *wdt;
  Waves *waves;
  static char *buf;
  int uidx, unit, unit0, n, mcnt, midx;

  if (argc < 3)
    return TCL_OK;
  wdt = (WDT *) atol (argv[1]);
  unit0 = atoi (argv[2]);
  mcnt = member_count[unit0];
  if (mcnt == 0)
    return TCL_OK;
  
  if (argc < 4 || (waves = (Waves *) atol (argv[3])) == 0)
    TCALLOC (waves, 1);

  TREALLOC (waves->y_ye, mcnt);
  TREALLOC (waves->rect, mcnt);
  waves->nsegments = 0;
  for (midx = 0; midx < mcnt; midx++) {
    unit = member[unit0][midx];
    for (uidx = 0; uidx < wdt->hdr->ccnt; uidx++)
      if (wdt->clus[uidx]->ccn == unit - 1)
	break;
    if (uidx == wdt->hdr->ccnt) {
      fprintf (stderr, "requested unit %d (zero based) is not in wdt file\n", unit - 1);
      fprintf (stderr, "units in file:");
      for (uidx = 0; uidx < wdt->hdr->ccnt; uidx++)
	fprintf (stderr, " %d", wdt->clus[uidx]->ccn);
      fprintf (stderr, "\n");
      return TCL_OK;
    }
    if (merge_map[unit] == unit) {
      unit == unit0 || DIE;
      memcpy (waves->v, wdt->clus[uidx]->raw, sizeof waves->v);
      for (n = 1; n < SNAPLEN; n++)
	if (waves->v[n] != waves->v[0])
	  break;
      if (n == SNAPLEN) {
	if (asprintf (&buf, "%d", unit) == -1) exit (1);
	Tcl_SetVar (interp, "unclassified_unit", buf, TCL_GLOBAL_ONLY | TCL_LIST_ELEMENT | TCL_APPEND_VALUE);
        free (buf);
      }
    }
    waves->y_ye[midx] = &wdt->clus[uidx]->y;
    {
      int tmp;
      char *p = (char *)&wdt->clus[uidx]->segcnt;
      memcpy (&tmp, p, sizeof tmp); /* some processors can't read unaligned int's */
      waves->nsegments += tmp;
    }
    waves->rect[midx] = &wdt->clus[uidx]->rect[0];
  }
  waves->member_count = mcnt;
  waves->unit0 = unit0;
  if (asprintf (&buf, "%ld", (long)waves) == -1) exit (1);
  Tcl_SetResult (interp, buf, (Tcl_FreeProc *)free);

  return TCL_OK;
}

static int
get_splot (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  SPL *spl;
  Splot *splot;
  static char *buf;
  int a, b, n, could_use;

  //  printf ("%s line %d: get_splot argc %d\n", __FILE__, __LINE__, argc);
  if (argc < 4)
    return TCL_OK;
  spl = (SPL *) atol (argv[1]);
  a = atoi (argv[2]) - 1;
  if ((b = atoi (argv[3])) > 0)
    b--;
  if (argc < 5 || (splot = (Splot *) atol (argv[4])) == 0)
    TCALLOC (splot, 1);

  could_use = -1;
  for (n = 0; n < spl->hdr->plane_count; n++) {
    if (spl->idx[n].a == a && spl->idx[n].b == b)
      break;
    if ((spl->idx[n].a == a || spl->idx[n].b == a) && b == -1)
      could_use = n;
  }

  if (0) {
    if (n == spl->hdr->plane_count)
      printf ("%s line %d: plane for %d and %d not found\n", __FILE__, __LINE__, a, b);
    else
      printf ("%s line %d: plane %d is %d and %d\n", __FILE__, __LINE__, n, a, b);
  }
  
  if (n == spl->hdr->plane_count && could_use >= 0)
    n = could_use;
  if (n == spl->hdr->plane_count)
    return TCL_OK;

  if (spl->plane[n] == 0) {
    (spl->plane[n] = malloc (spl->idx[n].count)) || DIE;
    fseek (spl->f, spl->idx[n].offset, SEEK_SET);
    fread (spl->plane[n], spl->idx[n].count, 1, spl->f) == 1 || DIE;
  }

  splot->a = a;
  splot->b = b;
  splot->plane = spl->plane[n];
  splot->narc = spl->narc[n];
  splot->used = spl->hdr->used;
  TREALLOC (splot->arclist, splot->used);
  TREALLOC (splot->arc, splot->narc);
  for (n = 0; n < splot->narc; n++) {
    splot->arc[n].width = 6;
    splot->arc[n].height = 6;
    splot->arc[n].angle1 = 0;
    splot->arc[n].angle2 = 368 * 64;
  }
  if (asprintf (&buf, "%ld", (long)splot) == -1) exit (1);
  Tcl_SetResult (interp, buf, (Tcl_FreeProc *)free);

  return TCL_OK;
}

static int
fromzero (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  SPL *spl;
  static char *buf;
  int a;
  float *d;

  if (argc < 3)
    return TCL_OK;
  if ((spl = (SPL *) atol (argv[1])) == 0)
    return TCL_OK;
  a = atoi (argv[2]);
  if (a >= spl->hdr->used - 1 || a < 0)
    return TCL_OK;
  d = &spl->hdr->d;
  if (asprintf (&buf, "%.1f", d[a]) == -1) exit (1);
  Tcl_SetResult (interp, buf, (Tcl_FreeProc *)free);
  
  return TCL_OK;
}

static int
get_t_hist(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Pos *pos;
  Hist *hist;
  static char *dounit;
  int bincnt, n, bin, max;
  double binstart, binend, binsize;
  double samples_per_ms = sampling_frequency / 1000;
  static char *buf;

  TREALLOC (dounit, id_count);
  TMEMSET (dounit, 0, id_count);
  if (argc < 6)
    return TCL_OK;
  pos = (Pos *) atol (argv[1]);
  binstart = atof (argv[2]);
  binend = atof (argv[3]);
  if (argc > 4)
    bincnt = atoi (argv[4]);
  else
    bincnt = 100;
  binsize = (binend - binstart) / bincnt;

  {
    int id, i;

    for (n = 6; n < argc; n++) {
      id =  atoi (argv[n]);
      (id > 0 && id < id_count) || DIE;
      for (i = 0; i < member_count[id]; i++)
	dounit[member[id][i]] = 1;
    }
  }
  
  hist = (Hist *) atol (argv[5]);
  if (hist) {
    if (bincnt == 0)
      bincnt = hist->bincnt;
    else {
      if (hist->bincnt != bincnt)
	TREALLOC (hist->bin, bincnt);
      hist->bincnt = bincnt;
    }
    memset (hist->bin, 0, bincnt * sizeof *hist->bin);
    if (binsize == 0)
      binsize = hist->binsize;
    else
      hist->binsize = binsize;
  }
  else {
    TCALLOC (hist, 1);
    TCALLOC (hist->bin, bincnt);
    hist->bincnt = bincnt;
    hist->binsize = binsize;
  }

  for (n = 0; n < pos->count; n++) {
    int id;
    double su = shift[pos->spike[n].id];
    double ms = (pos->spike[n].loc + su) / samples_per_ms;

    if ((id = pos->spike[n].id) >= 0 && id <= id_count && dounit[id]
	&& ms >= binstart && ms < binend) {
      bin = (ms - binstart) / binsize;
      if (bin < hist->bincnt)
	hist->bin[bin]++;
    }
  }

  for (max = n = 0; n < bincnt; n++)
    if (hist->bin[n] > max)
      max = hist->bin[n];
  hist->maxbin = max;
  if (asprintf (&buf, "%ld", (long)hist) == -1) exit (1);
  Tcl_SetResult (interp, buf, (Tcl_FreeProc *)free);

  return TCL_OK;
}

static int
get_members (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  static char buf[12];
  int n, i, id;

  for (n = 1; n < argc; n++) {
    id =  atoi (argv[n]);
    (id > 0 && id < id_count) || DIE;
    for (i = 0; i < member_count[id]; i++) {
      sprintf (buf, "%d", member[id][i]);
      Tcl_AppendElement (interp, buf);
    }
  }
  return TCL_OK;
}

static int
get_all (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  static char buf[12];
  int n;

  for (n = 0; n < uidx_count; n++) {
    sprintf (buf, "%d", unit_list[n]);
    Tcl_AppendElement (interp, buf);
  }
  return TCL_OK;
}

static int
get_hist(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Pos *pos;
  Hist *hist;
  static char *dounit;
  int bincnt, n, bin, max;
  double binsize, isi;
  double last = HUGE_VAL;
  double samples_per_ms = sampling_frequency / 1000;
  static char *buf;
  double binstart, binend;

  TREALLOC (dounit, id_count);
  TMEMSET (dounit, 0, id_count);
  if (argc < 6)
    return TCL_OK;
  pos = (Pos *) atol (argv[1]);
  binsize = atof (argv[2]);
  if (argc > 4)
    bincnt = atoi (argv[3]);
  else
    bincnt = 100;

  {
    int id, i;

    for (n = 5; n < argc; n++) {
      id =  atoi (argv[n]);
      (id > 0 && id < id_count) || DIE;
      for (i = 0; i < member_count[id]; i++)
	dounit[member[id][i]] = 1;
    }
  }
  
  hist = (Hist *) atol (argv[4]);
  if (hist) {
    if (bincnt == 0)
      bincnt = hist->bincnt;
    else {
      if (hist->bincnt != bincnt)
	TREALLOC (hist->bin, bincnt);
      hist->bincnt = bincnt;
    }
    memset (hist->bin, 0, bincnt * sizeof *hist->bin);
    if (binsize == 0)
      binsize = hist->binsize;
    else
      hist->binsize = binsize;
  }
  else {
    TCALLOC (hist, 1);
    TCALLOC (hist->bin, bincnt);
    hist->bincnt = bincnt;
    hist->binsize = binsize;
  }

  binstart = atoi (Tcl_GetVar (interp, "binstart", TCL_GLOBAL_ONLY));
  binend = atoi (Tcl_GetVar (interp, "binend", TCL_GLOBAL_ONLY));
  binstart *= samples_per_ms;
  binend *= samples_per_ms;

  int dblcmp (const void *va, const void *vb)
  {
    const double *a = (const double *) va;
    const double *b = (const double *) vb;

    if (*a == *b) return 0;
    if (*a < *b) return -1;
    return 1;
  }

  double *loc; int loccnt = 0;
  TMALLOC (loc, pos->count);
  for (n = 0; n < pos->count; n++) {
    int id;
    double su = shift[pos->spike[n].id];
    if ((id = pos->spike[n].id) >= 0 && id <= id_count && dounit[id]
	&& (pos->spike[n].loc + su) >= binstart
	&& (pos->spike[n].loc + su) <= binend) {
      loc[loccnt++] = pos->spike[n].loc + su;
    }
  }
  qsort (loc, loccnt, sizeof *loc, dblcmp);
  for (n = 0; n < loccnt; n++) {
    if (last != HUGE_VAL) {
      isi = (loc[n] - last) / samples_per_ms;
      bin = isi / binsize;
      if (bin < hist->bincnt)
        hist->bin[bin]++;
    }
    last = loc[n];
  }
  free (loc);

  for (max = n = 0; n < bincnt; n++)
    if (hist->bin[n] > max)
      max = hist->bin[n];
  hist->maxbin = max;
  if (asprintf (&buf, "%ld", (long)hist) == -1) exit (1);
  Tcl_SetResult (interp, buf, (Tcl_FreeProc *)free);
  
  return TCL_OK;
}

static int
open_pos(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  FILE *f;
  int count;
  Pos *pos;
  static char *buf;
  char *name;

  if (argc < 2)
    return TCL_OK;
  if ((f = fopen (name = change_filetype_trim (argv[1], ".pos", 0), "rb")) == 0) {
    fprintf (stderr, "error opening \"%s\": ", name);
    perror ("");
    return TCL_OK;
  }

  if (argc > 2)
    pos = (Pos *) atol (argv[2]);
  else
    pos = 0;

  count = filesize (f) / sizeof (Entry);
  filesize (f) == count *  sizeof (Entry) || DIE;
  if (pos)
    TREALLOC (pos->spike, count);
  else {
    TCALLOC (pos, 1);
    TMALLOC (pos->spike, count);
  }
  TREALLOC (pos->filename, strlen (name) + 1);
  strcpy (pos->filename, name);
  {
    FILE *g;
    struct stat s;
    int namelen;

    name = change_filetype_trim (argv[1], ".chan", 0);
    if (stat (name, &s) == 0) {
      TREALLOC (pos->chan_path, strlen (name) + 1);
      strcpy (pos->chan_path, name);
      goto read_pos;
    }
    name = change_filetype_trim (argv[1], ".path", 0);
    if ((g = fopen (name, "rb"))) {
      namelen = filesize (g);
      TREALLOC (pos->chan_path, namelen);
      fread (pos->chan_path, 1, namelen, g) == namelen || DIE;
      fclose (g);
      goto read_pos;
    }
    pos->chan_path = 0;
    goto read_pos;

    free (pos->spike);
    free (pos);
    return TCL_OK;
  }
 read_pos:
  fread (pos->spike, sizeof *pos->spike, count, f) == count || DIE;
  pos->count = count;
  if (asprintf (&buf, "%ld", (long)pos) == -1) exit (1);
  Tcl_SetResult (interp, buf, (Tcl_FreeProc *)free);
  return TCL_OK;
}

static int
open_wdt(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  FILE *f;
  int count, cidx;
  WDT *wdt;
  unsigned short *buf, *bufp, *endp;
  static char *ibuf;
  char *name;

  if (argc < 2)
    return TCL_OK;
  if ((f = fopen (name = change_filetype_trim (argv[1], ".wdt", 0), "rb")) == 0) {
    fprintf (stderr, "error opening \"%s\": ", name);
    perror ("");
    return TCL_OK;
  }
  
  if (argc > 2)
    wdt = (WDT *) atol (argv[2]);
  else
    wdt = 0;
  count = filesize (f) / sizeof *buf;
  if (wdt) {
    buf = (unsigned short *) wdt->hdr;
    TREALLOC (buf, count);
  }
  else {
    TCALLOC (wdt, 1);
    TMALLOC (buf, count);
  }
  wdt->hdr = (WDT_Hdr *)buf;
  fread (buf, sizeof *buf, count, f) == count || DIE;
  TREALLOC (wdt->clus, wdt->hdr->ccnt);
  endp = bufp = buf + count;
  for (cidx = wdt->hdr->ccnt - 1; cidx >= 0; cidx--) {
    bufp -= sizeof (unsigned) / sizeof *buf;
    memcpy (&count, bufp, sizeof count);
    //    count = *(unsigned *) bufp;
    bufp -= count;
    (bufp > buf && bufp < endp) || DIE;
    wdt->clus[cidx] = (Clus *)bufp;
    //    printf ("open_wdt cidx %d, ccn %d, count %d, segcnt %d\n", cidx, wdt->clus[cidx]->ccn, count, wdt->clus[cidx]->segcnt);
    //    printf ("cluster %d spikestart %d\n", cidx + 1, wdt->clus[cidx]->spikestart);
  }

  if (asprintf (&ibuf, "%ld", (long)wdt) == -1) exit (1);
  Tcl_SetResult (interp, ibuf, (Tcl_FreeProc *)free);
  return TCL_OK;
}

static int
open_spl (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  FILE *f;
  int n;
  SPL *spl;
  static char *ibuf;
  char *name;
  static int init_header_size = sizeof spl->hdr->plane_count + sizeof spl->hdr->used;
  int header_size;

  if (argc < 2)
    return TCL_OK;
  if ((f = fopen (name = change_filetype_trim (argv[1], ".spl", 0), "rb")) == 0) {
    fprintf (stderr, "error opening \"%s\": ", name);
    perror ("");
    return TCL_OK;
  }
  
  if (argc > 2)
    spl = (SPL *) atol (argv[2]);
  else
    spl = 0;
  
  if (spl) {
    for (n = 0; n < spl->hdr->plane_count; n++) {
      free (spl->plane[n]);
      spl->plane[n] = 0;
    }
    fclose (spl->f);
  }
  else {
    TCALLOC (spl, 1);
    (spl->hdr = malloc (init_header_size)) || DIE;
  }
    
  fread (spl->hdr, init_header_size, 1, f) == 1 || DIE;
  header_size = init_header_size + (spl->hdr->used - 1) * sizeof spl->hdr->d;
  TREALLOC (spl->hdr, header_size);
  rewind (f);
  fread (spl->hdr, header_size, 1, f) == 1 || DIE;

  TREALLOC (spl->plane, spl->hdr->plane_count);
  TMEMSET (spl->plane, 0, spl->hdr->plane_count);
  TREALLOC (spl->narc, spl->hdr->plane_count);
  TREALLOC (spl->idx, spl->hdr->plane_count);

  {
    struct {int count; short a, b;} tmp;
    int offset = filesize (f);
    fseek (f, -4, SEEK_END);
    fread (&tmp, sizeof tmp.count, 1, f) == 1 || DIE;
    for (n = spl->hdr->plane_count - 1; n >= 0; n--) {
      tmp.count > 0 || DIE;
      spl->narc[n] = tmp.count - 1 - 1 - 2 - 2 * spl->hdr->used; /* a (short), b (short), d (float), (ccn (short), 0x8000 (short)) for each used */
      spl->narc[n] / 2 * 2 == spl->narc[n] || DIE;
      spl->idx[n].count = tmp.count * 2 + 4;
      offset -= spl->idx[n].count;
      offset >= header_size || DIE;
      spl->idx[n].offset = offset;
      fseek (f, offset - 4, SEEK_SET);
      fread (&tmp, sizeof tmp, 1, f) == 1 || DIE;
      spl->idx[n].a = tmp.a;
      spl->idx[n].b = tmp.b;
    }
    offset == header_size || DIE;
  }
  spl->f = f;

  if (asprintf (&ibuf, "%ld", (long)spl) == -1) exit (1);
  Tcl_SetResult (interp, ibuf, (Tcl_FreeProc *)free);
  return TCL_OK;
}

static int
get_info(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Tk_Window tkwin = Tk_MainWindow(interp);
  Window w = Tk_WindowId(tkwin);
  Display *display = Tk_Display (tkwin);
  Window root_return, parent_return, *children_return;
  unsigned nchildren_return;
  XWindowAttributes window_attributes_return;
  int x_return, y_return;
  unsigned int width_return, height_return;
  unsigned int border_width_return;
  unsigned int depth_return;

  printf ("window id: %ld\n", w);

  XGetWindowAttributes(display, w, &window_attributes_return);
  printf ("width of %ld: %d %d\n", w, window_attributes_return.width, window_attributes_return.height);
  XGetGeometry(display, w, &root_return, &x_return, &y_return, &width_return, &height_return, &border_width_return, &depth_return);
  printf ("geometry of %ld: %d %d\n", w, width_return, height_return);

  XQueryTree(display, w, &root_return, &parent_return, &children_return, &nchildren_return);
  printf ("parent: %ld, root %ld\n", parent_return, root_return);

  XGetWindowAttributes(display, parent_return, &window_attributes_return);
  printf ("width of %ld: %d %d\n", parent_return, window_attributes_return.width, window_attributes_return.height);
  XGetGeometry(display, parent_return, &root_return, &x_return, &y_return, &width_return, &height_return, &border_width_return, &depth_return);
  printf ("geometry of %ld: %d %d\n", parent_return, width_return, height_return);

  XQueryTree(display, parent_return, &root_return, &parent_return, &children_return, &nchildren_return);
  printf ("parent: %ld, root %ld\n", parent_return, root_return);

  XGetWindowAttributes(display, parent_return, &window_attributes_return);
  printf ("width of %ld: %d %d\n", parent_return, window_attributes_return.width, window_attributes_return.height);
  XGetGeometry(display, parent_return, &root_return, &x_return, &y_return, &width_return, &height_return, &border_width_return, &depth_return);
  printf ("geometry of %ld: %d %d\n", parent_return, width_return, height_return);

  XQueryTree(display, parent_return, &root_return, &parent_return, &children_return, &nchildren_return);
  printf ("parent: %ld, root %ld\n", parent_return, root_return);
  
  return TCL_OK;
}

static int
merge (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  int n, id, id0, midx;
  Pos *pos;

  (pos = (Pos *) atol (argv[1])) || DIE;
  argc > 1 || DIE;
  new_mrg (interp, pos);
  if (argc == 2) {
    for (n = 0; n < id_count; n++)
      if (merge_map[n] != -1)
	merge_map[n] = n;
  }
  else {
    id0 = atoi (argv[2]);
    (id0 > 0 && id0 < id_count) || DIE;
    member_count[id0] || DIE;
    if (argc == 3) {
      for (n = 0; n < member_count[id0]; n++) {
	id = member[id0][n];
	merge_map[id] = id;
      }
    }
    else { /* argc > 3 */
      for (n = 2; n < argc; n++) {
	id = atoi (argv[n]);
	(id > 0 && id < id_count) || DIE;
	for (midx = 0; midx < member_count[id]; midx++)
	  merge_map[member[id][midx]] = id0;
      }
    }
  }
  write_mrg (interp, pos);
  match_map (interp, pos);

  return TCL_OK;
}

static int
delete (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  int n, id, id0, midx;
  Pos *pos;

  (pos = (Pos *) atol (argv[1])) || DIE;
  argc > 1 || DIE;
  new_mrg (interp, pos);
  if (argc == 2) {
    for (n = 0; n < id_count; n++)
      if (merge_map[n] == -1)
	merge_map[n] = n;
  }
  else {
    for (n = 2; n < argc; n++) {
      id0 = atoi (argv[n]);
      (id0 > 0 && id0 < id_count) || DIE;
      for (midx = 0; midx < member_count[id0]; midx++) {
	id = member[id0][midx];
	merge_map[id] = -1;
      }
    }
  }
  write_mrg (interp, pos);
  match_map (interp, pos);

  return TCL_OK;
}

static int
get_doublet (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Pos *pos;
  static char *dounit;
  int n;
  double isi, last;
  double samples_per_ms = sampling_frequency / 1000;

  if (argc < 3)
    return TCL_OK;
  TREALLOC (dounit, id_count);
  TMEMSET (dounit, 0, id_count);
  pos = (Pos *) atol (argv[1]);

  {
    int id, i;

    for (n = 2; n < argc; n++) {
      id =  atoi (argv[n]);
      (id > 0 && id < id_count) || DIE;
      for (i = 0; i < member_count[id]; i++)
	dounit[member[id][i]] = 1;
    }
  }

  double su = shift[pos->spike[pos->doublet_start].id];
  last = pos->spike[pos->doublet_start].loc + su;
  for (n = pos->doublet_start + 1; n != pos->doublet_start; n = (n + 1) % pos->count) {
    int id;
    double su = shift[pos->spike[n].id];
    if ((id = pos->spike[n].id) >= 0 && id <= id_count && dounit[id]) {
      isi = ((pos->spike[n].loc + su) - last) / samples_per_ms;
      if (isi >= 0 && isi < 3)
	break;
      last = pos->spike[n].loc + su;
    }
  }
  if (n == pos->doublet_start)
    return TCL_OK;
  pos->doublet_start = n;

  {
    static char *cmd;
    char *filename = change_filetype_trim (pos->filename, ".chan", 0);
    char *fmt = "waveform.tcl %s %.0f %.0f";

    double su = shift[pos->spike[n].id];
    if (asprintf (&cmd, fmt, filename, last - 20, pos->spike[n].loc + su + 44) == -1) exit (1);
    Tcl_SetResult (interp, cmd, (Tcl_FreeProc *)free);
  }
  
  return TCL_OK;
}

static int
pos_retype (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Pos *pos;
  char *fname;

  if (argc < 3)
    return TCL_OK;
  pos = (Pos *) atol (argv[1]);
  fname = change_filetype_trim (pos->filename, argv[2], 0);
  Tcl_AppendElement (interp, fname);
  return TCL_OK;
}

static int
merge_version (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  FILE *f;
  Pos *pos;
  char *mrg, *mrgv;

  argc == 3 || DIE;
  (pos = (Pos *) atol (argv[1])) || DIE;

  new_mrg (interp, pos);
  if (argv[2][0] == '0') {
    int n;
    for (n = 0; n < id_count; n++)
      merge_map[n] = n;
    write_mrg (interp, pos);
    return TCL_OK;
  }
  mrg = change_filetype_trim (pos->filename, ".mrg", 0);
  if (asprintf (&mrgv, "%s.~%s~", mrg, argv[2]) == -1) exit (1);

  (f = fopen (mrgv, "rb")) || DIE;
  free (mrgv);
  fread (merge_map, sizeof *merge_map, id_count, f) == id_count || DIE;
  fclose (f);

  write_mrg (interp, pos);

  return TCL_OK;
}

static int
ymag (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  argc == 2 || DIE;
  yscale_mult = atof (argv[1]);
  return TCL_OK;
}

static int
set_shiftpxl (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  argc == 2 || DIE;
  shiftpxl = atoi (argv[1]);
  return TCL_OK;
}

int
Tcl_AppInit(Tcl_Interp *interp)
{
  if (Tcl_Init (interp) == TCL_ERROR)
    return TCL_ERROR;

  if (Tk_Init (interp) == TCL_ERROR)
    return TCL_ERROR;


  Tcl_CreateCommand (interp, "square", SquareCmd, (ClientData) Tk_MainWindow(interp), 0);
  Tcl_CreateCommand (interp, "putint", PutIntCmd, 0, 0);
  Tcl_CreateCommand (interp, "write_edt", write_edt, 0, 0);
  Tcl_CreateCommand (interp, "open_pos", open_pos, 0, 0);
  Tcl_CreateCommand (interp, "open_wdt", open_wdt, 0, 0);
  Tcl_CreateCommand (interp, "open_spl", open_spl, 0, 0);
  Tcl_CreateCommand (interp, "get_units", get_units, 0, 0);
  Tcl_CreateCommand (interp, "get_t_hist", get_t_hist, 0, 0);
  Tcl_CreateCommand (interp, "get_hist", get_hist, 0, 0);
  Tcl_CreateCommand (interp, "get_waves", get_waves, 0, 0);
  Tcl_CreateCommand (interp, "get_splot", get_splot, 0, 0);
  Tcl_CreateCommand (interp, "get_maxbin", get_maxbin, 0, 0);
  Tcl_CreateCommand (interp, "get_last_time", get_last_time, 0, 0);
  Tcl_CreateCommand (interp, "get_color", get_color, 0, 0);
  Tcl_CreateCommand (interp, "fromzero", fromzero, 0, 0);
  Tcl_CreateCommand (interp, "get_distance", get_distance, 0, 0);
  Tcl_CreateCommand (interp, "get_info", get_info, 0, 0);
  Tcl_CreateCommand (interp, "get_mwaves", get_mwaves, 0, 0);
  Tcl_CreateCommand (interp, "waitpid", wait_pid, 0, 0);
  Tcl_CreateCommand (interp, "merge", merge, 0, 0);
  Tcl_CreateCommand (interp, "delete", delete, 0, 0);
  Tcl_CreateCommand (interp, "get_members", get_members, 0, 0);
  Tcl_CreateCommand (interp, "get_all", get_all, 0, 0);
  Tcl_CreateCommand (interp, "get_doublet", get_doublet, 0, 0);
  Tcl_CreateCommand (interp, "pos_retype", pos_retype, 0, 0);
  Tcl_CreateCommand (interp, "merge_version", merge_version, 0, 0);
  Tcl_CreateCommand (interp, "ymag", ymag, 0, 0);
  Tcl_CreateCommand (interp, "set_shiftpxl", set_shiftpxl, 0, 0);
  Tcl_CreateCommand (interp, "write_shift", write_shift, 0, 0);
  Tcl_CreateCommand (interp, "recolor", recolor, 0, 0);

  return TCL_OK;
}

int main(int argc, char **argv)
{
  Tk_Main(argc, argv, Tcl_AppInit);
  return 0;
}


