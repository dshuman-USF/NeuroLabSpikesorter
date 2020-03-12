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


/* window.c */


#include "nde.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <string.h>

static Display *display;
static int screen;
static Window win;
static unsigned int width, height;
static XFontStruct *font_info;
static GC gc;
static int mapped, changed = 1;
static double yscale = 32768;
static double xclick, yclick;
static int nowait;

static
struct {
  /*@observer@*/  char *color_name;
  unsigned long color;
  float v[SNAPLEN];
} *snaps;

static int snaps_count, snaps_alloc;
static char *text;
static int textlen;

void
get_click (double *xp, double *yp)
{
  *xp = xclick;
  *yp = yclick;
  return;
}

void
draw_snap (float *v, char *color)
{
  int snap_idx = snaps_count++;

  if (snaps_count > snaps_alloc)
    TREALLOC (snaps, snaps_alloc = snaps_count);
  memmove (snaps[snap_idx].v, v, sizeof snaps[snap_idx].v);
  snaps[snap_idx].color_name = color;
  changed = 1;
}

void
draw_text (char *t)
{
  textlen =  (int)strlen (t);
  TREALLOC (text, textlen + 1);
  strcpy (text, t);
}

void
draw_peak (float *v, char *color)
{
  int cnt = 2 * WRADIUS0 + 1, n;
  float *w;
  int snap_idx = snaps_count++;

  if (snaps_count > snaps_alloc)
    TREALLOC (snaps, snaps_alloc = snaps_count);
  memset (snaps[snap_idx].v, 0, sizeof snaps[snap_idx].v);
  memmove (snaps[snap_idx].v + PRESAMPLES - WRADIUS0, v - WRADIUS0, (2 * WRADIUS0 + 1) * sizeof (float));
  w = snaps[snap_idx].v + PRESAMPLES + 2 * WRADIUS0;
  memcpy (w, v - WRADIUS0, cnt * sizeof (float));
  whiten_samples (w, cnt);
  for (n = 0; n < cnt; n++)
    w[n] *= 3200;
  snaps[snap_idx].color_name = color;
  changed = 1;
}

static void
show_text ()
{
  static int    font_height;
  static int    msg_x, msg_y;

  XSetForeground(display, gc, BlackPixel (display, screen));
  XSetBackground(display, gc, WhitePixel (display, screen));
  XTextWidth (font_info, text, textlen);
  font_height = font_info->ascent + font_info->descent;
  msg_x  = font_height;
  msg_y  = 2 * font_height;
  XDrawImageString (display, win, gc, msg_x, msg_y, text, textlen);
}

static void
show_snaps (void)
{
  int n, i;
  XPoint p[SNAPLEN];
  double ymin, ymax, ydif, x, y;
  XColor xcolor;
  static /*@observer@*/ char *gray_name = "light gray";
  static unsigned long gray;

  if (gray_name) {
    XAllocNamedColor (display, DefaultColormap (display, screen), gray_name, &xcolor, &xcolor);
    gray_name = 0;
    gray = xcolor.pixel;
  }
  ymin = -fabs (yscale);
  ymax = fabs (yscale);
  ydif = ymax - ymin;
  XSetForeground(display, gc, gray);
  x = (double)PRESAMPLES / (SNAPLEN - 1) * (width - 1);
  XDrawLine (display, win, gc, (int)x, 0, (int)x, (int)(height - 1));

  for (i = 0; i < SNAPLEN; i++) {
    double f = i % 10 == 0 ? 0.0 : 0.0;
    x = (double)i / (SNAPLEN - 1) * (width - 1);
    if (i % 10 != 0)
      XDrawLine (display, win, gc, (int)x, (int)(f * height), (int)x, (int)((1 - f) * height));
  }

  XSetForeground(display, gc, BlackPixel(display, screen));
  for (i = 0; i < SNAPLEN; i++) {
    double f = i % 10 == 0 ? 0.0 : 0.0;
    x = (double)i / (SNAPLEN - 1) * (width - 1);
    if (i % 10 == 0)
      XDrawLine (display, win, gc, (int)x, (int)(f * height), (int)x, (int)((1 - f) * height));
  }

  for (n = 0; n < snaps_count; n++) {
    if (snaps[n].color_name) {
      XAllocNamedColor (display, DefaultColormap (display, screen), snaps[n].color_name, &xcolor, &xcolor);
      snaps[n].color_name = 0;
      snaps[n].color = xcolor.pixel;
    }
    for (i = 0; i < SNAPLEN; i++) {
      p[i].x = (short)((double)i / (SNAPLEN - 1) * (width - 1));
      y =  snaps[n].v[i];
      if (y > ymax) y = ymax;
      if (y < ymin) y = ymin;
      p[i].y = (short)((ymax - y) / ydif * (height - 1));
    }
    XSetForeground(display, gc, snaps[n].color);
    XDrawLines (display, win, gc, p, SNAPLEN, CoordModeOrigin);
  }
}

static void
show_all (void)
{
  show_snaps ();
  show_text ();
  changed = 0;
}

void clear (void)
{
  snaps_count = 0;
  changed = 1;
}

void
wclose (void)
{
  XUnloadFont(display, font_info->fid);
  XFreeGC(display, gc);
  XCloseDisplay(display);
  mapped = 0;
  display = 0;
}

#define BUFSZ 8

void window (void)
{
  int          x, y;
  unsigned int border_width = 0;
  XGCValues     values;

  if (display)
    return;
  (display = XOpenDisplay (0)) || DIE;

  x = y = 0;
  width  = 640;
  height = 480;

  win = XCreateSimpleWindow(display, RootWindow(display, screen),
			    x, y, width, height, border_width,
			    BlackPixel(display, screen),
			    WhitePixel(display, screen));

  XSelectInput(display, win, ExposureMask | KeyPressMask | ButtonPressMask | StructureNotifyMask);

  (font_info = XLoadQueryFont(display, "9x15")) || DIE;

  gc = XCreateGC (display, win, 0, &values);

  XSetFont(display, gc, font_info->fid);
  XSetForeground(display, gc, BlackPixel(display, screen));

  XMapWindow(display, win);
}

static int
handle_button (XButtonEvent *ev)
{
  double ymin, ymax, ydif;

  ymin = -fabs (yscale);
  ymax = fabs (yscale);
  ydif = ymax - ymin;

  xclick = (double) ev->x / (width - 1) * (SNAPLEN - 1);
  yclick = ymax - (double) ev->y / (height - 1) * ydif;;

  return (int)ev->button;
}

static int
event_loop ()
{
  XEvent report;
  static char buf[BUFSZ];
  static int escape;

  while (1) {
    if (nowait && QLength (display) == 0) {
      XFlush (display);
      return 0;
    }
    XNextEvent(display, (XEvent *)&report.xany);

    switch (report.type) {
    case Expose:
      if ( report.xexpose.count != 0 )
	break;
      mapped = 1;
      show_all ();
      break;
    case ConfigureNotify:
      width  = (unsigned)report.xconfigure.width;
      height = (unsigned)report.xconfigure.height;
      XClearWindow (display, win);
      show_all ();
      break;
    case MappingNotify:
      if (report.xmapping.request == MappingKeyboard || report.xmapping.request == MappingModifier)
	XRefreshKeyboardMapping (&report.xmapping);
      break;
    case ButtonPress:
      return handle_button (&report.xbutton);
    case KeyPress:
      XLookupString (&report.xkey, buf, BUFSZ, 0, 0);
      if (buf[0] == 27)
	escape = 1;
      else if (buf[0]) {
	if (escape) {
	  escape = 0;
	  return buf[0] + 128;
	}
	else return buf[0];
      }
    }
  }
}

int show (void)
{
  if (mapped && changed) {
    XClearWindow (display, win);
    show_all ();
  }
  return event_loop ();
}

int show_nowait (void)
{
  int v;

  nowait = 1;
  v = show ();
  nowait = 0;
  return v;
}

int
xor_snap (float *v, int wait)
{
  int i;
  XPoint p[SNAPLEN];
  double yscale = 32768;
  double ymin, ymax, ydif, y;
  XColor xcolor;
  static /*@observer@*/ char *color_name = "cyan";
  static unsigned long color = 0;

  if (color_name) {
    XAllocNamedColor (display, DefaultColormap (display, screen), color_name, &xcolor, &xcolor);
    color_name = 0;
    color = xcolor.pixel;
  }

  if (yscale != 0) {
    ymin = -fabs (yscale);
    ymax = fabs (yscale);
  }
  ydif = ymax - ymin;

  for (i = 0; i < SNAPLEN; i++) {
    p[i].x = (short)((double)i / (SNAPLEN - 1) * (width - 1));
    y =  v[i];
    if (y > ymax) y = ymax;
    if (y < ymin) y = ymin;
    p[i].y = (short)((ymax - y) / ydif * (height - 1));
  }
  XSetForeground(display, gc, color);
  XSetFunction(display, gc, GXxor);
  XDrawLines (display, win, gc, p, SNAPLEN, CoordModeOrigin);
  XSetFunction(display, gc, GXcopy);
  show_text ();
  if (wait)
    return event_loop ();
  else
    return 0;
}
