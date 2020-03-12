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

/* waveform.c */

#define _GNU_SOURCE 1

#include <config.h>
#include "tmalloc.h"
#include <math.h>
#include <tk.h>
#include <tcl.h>
#include <string.h>
#include <limits.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <ctype.h>
#include <X11/Xutil.h>

#ifdef HAVE_AO_PLAY
#include <ao/ao.h>
#endif

#define SF 25000

typedef struct {
  FILE *f;
  XPoint *point;
  int point_count, point_alloc;
  size_t sample_count, samples_left, samples_per_width;
  double xscale, yscale;
  int y0;
  long long starting_pixel;
  size_t next_sample;
  int in, out, max, min, last_x;
} Chan;

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
  Tk_Window tkwin;		/* Window that embodies the square.  NULL
				 * means window has been deleted but
				 * widget record hasn't been cleaned up yet. */
  Display *display;		/* X's token for the window's display. */
  Tcl_Interp *interp;		/* Interpreter associated with widget. */
  Tcl_Command widgetCmd;	/* Token for square's widget command. */
  int run;
  double time, xscale, yscale;
  Tcl_TimerToken ttoken, ttoken2;
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
} Square;

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
  {TK_CONFIG_DOUBLE, "-xscale", "xscale", "Xscale", "0", Tk_Offset(Square, xscale), 0},
  {TK_CONFIG_DOUBLE, "-yscale", "yscale", "Yscale", "0", Tk_Offset(Square, yscale), 0},
  {TK_CONFIG_INT, "-ttoken", "ttoken", "Ttoken", "0", Tk_Offset(Square, ttoken), 0},
  {TK_CONFIG_INT, "-ttoken2", "ttoken2", "Ttoken2", "0", Tk_Offset(Square, ttoken2), 0},
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
static void SquareEventProc _ANSI_ARGS_((ClientData clientData, XEvent *eventPtr));
static int SquareWidgetCmd _ANSI_ARGS_((ClientData clientData, Tcl_Interp *, int argc, const char **argv));
static void WFDisplay (ClientData clientData) ;
static void WFConfigure (ClientData clientData) ;
void nullproc (ClientData x) {}

void (*paint[])(ClientData) = {0,0,0,0,WFDisplay};
void (*resize[])(ClientData) = {0,0,0,0,WFConfigure};


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
  squarePtr->updatePending = 0;
  squarePtr->color = 0;
  squarePtr->ttoken = 0;
  squarePtr->ttoken2 = 0;
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

static void
WFDisplay (ClientData clientData)
{
  Square *squarePtr = (Square *) clientData;
  Tk_Window tkwin = squarePtr->tkwin;
  Pixmap pm = None;
  Drawable d;
  Chan *chan;
  static GC gc_red;

  squarePtr->updatePending = 0;
  if (!Tk_IsMapped(tkwin))
    return;

  if (squarePtr->doubleBuffer) {
    pm = Tk_GetPixmap(Tk_Display(tkwin), Tk_WindowId(tkwin), Tk_Width(tkwin), Tk_Height(tkwin), DefaultDepthOfScreen(Tk_Screen(tkwin)));
    d = pm;
  }
  else d = Tk_WindowId(tkwin);

  if (gc_red == 0) gc_red = Tk_GCForColor (Tk_GetColor (squarePtr->interp, tkwin, "red"), d);

  Tk_Fill3DRectangle(tkwin, d, squarePtr->bgBorder, 0, 0, Tk_Width(tkwin), Tk_Height(tkwin), squarePtr->borderWidth, squarePtr->relief);

  if ((chan = (Chan *)squarePtr->data)) {
    int yhi = chan->y0 - floor ((32767 * chan->yscale) + .5);
    int ylo = chan->y0 - floor ((-32768 * chan->yscale) + .5);
    XDrawLine (Tk_Display(tkwin), d, gc_red, 0, chan->y0, Tk_Width(tkwin)-1, chan->y0);
    XDrawLine (Tk_Display(tkwin), d, gc_red, 0, yhi, Tk_Width(tkwin)-1, yhi);
    XDrawLine (Tk_Display(tkwin), d, gc_red, 0, ylo, Tk_Width(tkwin)-1, ylo);
    XDrawLines (Tk_Display(tkwin), d, squarePtr->gc, chan->point, chan->point_count, CoordModeOrigin);
  }

  if (squarePtr->doubleBuffer) {
    XCopyArea(Tk_Display(tkwin), pm, Tk_WindowId(tkwin), squarePtr->gc, 0, 0, (unsigned) Tk_Width(tkwin), (unsigned) Tk_Height(tkwin), 0, 0);
    Tk_FreePixmap(Tk_Display(tkwin), pm);
  }
}

static char *
change_filetype_trim (const char *file_name, char *suffix, int trim)
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

#define BYTES_PER_SAMPLE 2

static long
filesize (FILE *f)
{
  long bytes;

  fseek (f, 0, SEEK_END);
  bytes = ftell (f);
  rewind (f);
  return bytes;
}

static int
open_chan(ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  FILE *f;
  Chan *chan;
  static char *buf;

  if (argc < 2)
    return TCL_OK;
  (f = fopen (change_filetype_trim (argv[1], ".chan", 0), "r")) || DIE;

  if (argc > 2)
    chan = (Chan *) atol (argv[2]);
  else
    chan = 0;
  if (!chan)
    TCALLOC (chan, 1);

  chan->f = f;
  chan->sample_count = filesize (f) / BYTES_PER_SAMPLE;

  if (asprintf (&buf, "%ld", (long)chan) == -1) exit (1);

  Tcl_SetResult (interp, buf, (Tcl_FreeProc *)free);
  return TCL_OK;
}

static inline long long
sample_to_pixel (Chan *chan, size_t sample)
{
  return floor (sample * chan->xscale + .5);
}

static inline size_t
pixel_to_prev_sample (Chan *chan, long long pixel)
{
  size_t sample;

  if (pixel < .5)
    return 0;
  sample = floor ((pixel - .5) / chan->xscale);
  
  while (sample > 0 && sample_to_pixel (chan, sample) == pixel)
    sample--;
  return sample;
}

static inline size_t
pixel_to_next_sample (Chan *chan, long long pixel)
{
  size_t sample = ceil ((pixel + .5) / chan->xscale);
  while (sample_to_pixel (chan, sample) == pixel)
    sample++;
  return sample;
}

#define BUFCNT 131072
#define ADD_POINT(xval,yval) (void)(pidx < chan->point_alloc || DIE), point[pidx].x = xval, point[pidx++].y = yval
#undef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

static void
get_next_chunk (ClientData clientData)
{
  Square *squarePtr = (Square *) clientData;
  Tk_Window tkwin = squarePtr->tkwin;
  Chan *chan = (Chan *)squarePtr->data;
  static short buf[BUFCNT];
  size_t readcnt;
  long long x;
  int y, n;
  XPoint *point = chan->point;
  long long starting_pixel = chan->starting_pixel;
  int last_x = chan->last_x;
  int in = chan->in;
  int out = chan->out;
  int max = chan->max;
  int min = chan->min;
  int pidx = chan->point_count;
  size_t next_sample = chan->next_sample;
  int point_pending = 0;
  unsigned short v = 1;

  Tcl_DeleteTimerHandler (squarePtr->ttoken);
  squarePtr->ttoken = 0;

  if (chan->samples_per_width > BUFCNT * 2 && atoi (Tcl_GetVar (squarePtr->interp, "nopaint", TCL_GLOBAL_ONLY))) {
    squarePtr->ttoken = Tcl_CreateTimerHandler(100, get_next_chunk, clientData);
    return;
  }

  readcnt = MIN (BUFCNT, chan->samples_left);
  fread (buf, sizeof *buf, readcnt, chan->f) == readcnt || DIE;
  if (*(char *)&v == 0) {
    /* big endian */
    size_t scnt = sizeof *buf * readcnt / 2;
    unsigned short *sbuf = (unsigned short *) buf;
    int n;

    for (n = 0; n < scnt; n++)
      sbuf[n] = ((sbuf[n] >> 8) & 0xff) | ((sbuf[n] & 0xff) << 8);
  }
  
  for (n = 0; n < readcnt; n++) {
    x = sample_to_pixel (chan, next_sample++);
    if (pidx + 1 > (x + 1) * 4 + 1) {
      printf ("pidx %d, x %lld, next_sample: %u\n", pidx, x, (unsigned)(next_sample - 1));
      exit (1);
    }
    x -= starting_pixel;

    y = chan->y0 - floor ((buf[n] * chan->yscale) + .5);
    if (x == last_x) {
      if      (y > max) max = y;
      else if (y < min) min = y;
      out = y;
      point_pending = 1;
    }
    else {
      if (max > in) ADD_POINT (last_x, max);
      if (min < in) ADD_POINT (last_x, min);
      if (out != min && out != max) ADD_POINT (last_x, out);
      ADD_POINT (x, y);
      in = out = min = max = y;
      point_pending = 0;
    }
    last_x = x;
  }
  chan->samples_left -= readcnt;
  if (chan->samples_left == 0) {
    if (point_pending) {
      if (max > in) ADD_POINT (last_x, max);
      if (min < in) ADD_POINT (last_x, min);
    }
  }
  else
    squarePtr->ttoken = Tcl_CreateTimerHandler(0, get_next_chunk, clientData);

  if (chan->samples_per_width > BUFCNT * 2) {
    int first = chan->point_count - 1;
    int count = pidx - first;

    if (first < 0) first = 0, count--;
    XDrawLines (Tk_Display(tkwin), Tk_WindowId(tkwin), squarePtr->gc, chan->point + first, count, CoordModeOrigin);
  }

  chan->last_x = last_x;
  chan->in = in;
  chan->out = out;
  chan->max = max;
  chan->min = min;
  chan->point_count = pidx;
  chan->next_sample = next_sample;
}

static inline long long
get_time (void)
{
  struct timeval tv;
  gettimeofday (&tv, 0);  
  return (long long) tv.tv_sec * 1000000 + tv.tv_usec;
}

#define BUFTIME 100000LL
#define BUFSAMP (SF * BUFTIME / 1000000)

#ifdef HAVE_AO_PLAY
static void
update_sound (ClientData clientData)
{
  Square *squarePtr = (Square *) clientData;
  Chan *chan = (Chan *)squarePtr->data;
  Tk_Window tkwin = squarePtr->tkwin;
  static ao_device *device;
  static ao_sample_format format;
  static int default_driver;
  int req_size, last_req = 0;
  static GC invert_gc;

  Tcl_DeleteTimerHandler (squarePtr->ttoken2);
  squarePtr->ttoken2 = 0;
  if (squarePtr->ttoken) {
    squarePtr->ttoken2 = Tcl_CreateTimerHandler (100, update_sound, clientData);
    return;
  }

  if (invert_gc == 0 && chan) {
    XGCValues values;
    values.function = GXinvert;
    invert_gc = XCreateGC (squarePtr->display, Tk_WindowId (tkwin), GCFunction, &values);
    format.bits = 16;
    format.channels = 1;
    format.rate = SF;
    format.byte_format = AO_FMT_LITTLE;
    ao_initialize();
    default_driver = ao_default_driver_id();
  }

  if (chan) {
    static long long start_time, sound_start_time;
    static size_t first_sample, next_sample;
    static short buffer[BUFSAMP];
    static int last_pixel = -1, now_pixel;
    size_t now_sample, end_sample, sound_now_sample;
    long long now;
    
    if (squarePtr->run == 1) {
      (device = ao_open_live(default_driver, &format, NULL /* no options */)) || DIE;
      next_sample = first_sample = pixel_to_prev_sample (chan, chan->starting_pixel);
      chan->f || DIE;
      fseek (chan->f, first_sample * BYTES_PER_SAMPLE, SEEK_SET);
      start_time = get_time ();
      sound_start_time = start_time;
      squarePtr->run = 2;
    }
    now = get_time ();

    now_sample = first_sample + (now - start_time) * SF / 1000000;
    if (now > sound_start_time)
      sound_now_sample = first_sample + (now - sound_start_time) * SF / 1000000;
    else
      sound_now_sample = first_sample;
    if (squarePtr->run == 0 || now_sample >= first_sample + chan->samples_per_width) {
      if (device) {
	ao_close(device);
	device = 0;
      }
      first_sample = next_sample = 0;
      last_pixel = -1; now_pixel = 0;
      start_time = sound_start_time = 0;
      squarePtr->run = 0;
      return;
    }
    
    now_pixel = (sound_now_sample - first_sample) * chan->xscale;
    if (now_pixel >= 0 && now_pixel > last_pixel) {
      if (last_pixel >= 0) {
	XDrawLine (Tk_Display(tkwin), Tk_WindowId (tkwin), invert_gc, last_pixel, 0, last_pixel, Tk_Height (tkwin));
      }
      XDrawLine (Tk_Display(tkwin), Tk_WindowId (tkwin), invert_gc, now_pixel, 0, now_pixel, Tk_Height (tkwin));
      last_pixel = now_pixel;
    }
    end_sample = now_sample + BUFSAMP;
    if (end_sample > chan->sample_count) end_sample = chan->sample_count, last_req = 1;
    if (end_sample > first_sample + chan->samples_per_width)
      end_sample = first_sample + chan->samples_per_width, last_req = 1;
    req_size = end_sample - next_sample;
    if (req_size > BUFSAMP) req_size = BUFSAMP;
    if (device && (req_size > BUFSAMP / 10 || last_req)) {
      fread (buffer, sizeof *buffer, req_size, chan->f) == req_size || DIE;
      ao_play (device, (char *)buffer, req_size * BYTES_PER_SAMPLE);
      next_sample += req_size;
      if (last_req) {
	ao_close(device);
	device = 0;
      }
    }
    squarePtr->ttoken2 = Tcl_CreateTimerHandler (0, update_sound, clientData);
  }
}
#endif

static void
WFConfigure (ClientData clientData)
{
  Square *squarePtr = (Square *) clientData;
  Tk_Window tkwin = squarePtr->tkwin;
  Chan *chan;
  int width = Tk_Width(tkwin);
  int height = Tk_Height(tkwin);

  if (width == 1)
    return;

  chan = (Chan *)squarePtr->data;

  if (chan) {
    size_t first_sample, last_sample, sample_count;

    chan->xscale = squarePtr->xscale;
    chan->yscale = squarePtr->yscale;
    if (!chan->yscale) chan->yscale = .005;
    if (!chan->xscale) chan->xscale = 1;
    
    chan->starting_pixel = sample_to_pixel (chan, (int)squarePtr->time);
    first_sample = pixel_to_prev_sample (chan, chan->starting_pixel);
    last_sample = pixel_to_next_sample (chan, chan->starting_pixel + width - 1);
    chan->f || DIE;
    sample_count = chan->sample_count;
    if (last_sample >= sample_count) last_sample = sample_count - 1;
    if (first_sample >= sample_count) first_sample = sample_count;
    fseek (chan->f, first_sample * BYTES_PER_SAMPLE, SEEK_SET);
    TREALLOC (chan->point, chan->point_alloc = width * 4 + 2);
    chan->point_count = 0;
    chan->samples_per_width = chan->samples_left = last_sample - first_sample + 1;
    chan->y0 = height / 2;
    chan->next_sample = first_sample;
    chan->last_x = sample_to_pixel (chan, first_sample) - chan->starting_pixel - 1;
    chan->in = chan->out = chan->min = chan->max = 0;

    if (chan->samples_per_width > BUFCNT * 2 && atoi (Tcl_GetVar (squarePtr->interp, "nopaint", TCL_GLOBAL_ONLY)) == 0)
      Tk_Fill3DRectangle(tkwin, Tk_WindowId(tkwin), squarePtr->bgBorder, 0, 0,
			 Tk_Width(tkwin), Tk_Height(tkwin), squarePtr->borderWidth, squarePtr->relief);

    if (!squarePtr->ttoken)
      get_next_chunk (clientData);

#ifdef HAVE_AO_PLAY
    if (!squarePtr->ttoken2 && squarePtr->run == 1)
      update_sound (clientData);
#endif

  }
}

static int
get_sample_count (ClientData clientData, Tcl_Interp *interp, int argc, const char **argv)
{
  Chan *chan;
  static char*buf;

  if (argc < 2)
    return TCL_OK;

  chan = (Chan *) atol (argv[1]);
  if (asprintf (&buf, "%u", (unsigned)chan->sample_count) == -1) exit (1);

  Tcl_SetResult (interp, buf, (Tcl_FreeProc *)free);
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
  Tcl_CreateCommand (interp, "open_chan", open_chan, 0, 0);
  Tcl_CreateCommand (interp, "get_sample_count", get_sample_count, 0, 0);

  return TCL_OK;
}

int main(int argc, char **argv)
{
  Tk_Main(argc, argv, Tcl_AppInit);
  return 0;
}
