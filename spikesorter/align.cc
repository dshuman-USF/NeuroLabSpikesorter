//g++ -Wall -g -O2 -std=gnu++11 `pkg-config --cflags gtk+-3.0` -o align align.cc `pkg-config --libs gtk+-3.0 gsl`
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <errno.h>
#include <error.h>
#include <gtk/gtk.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;

class Spline
{
  gsl_interp_accel *acc;
  gsl_spline *spline;
public:
  Spline():acc(0), spline(0) {};
  Spline(array<short, 65> &sv);
  ~Spline();
  double val (double x) { return gsl_spline_eval (spline, x, acc);}
};

Spline::Spline (array<short, 65> &sv)
{
  vector<double> xv, yv;
  for (size_t i = 0; i < sv.size(); ++i) {
    xv.push_back (i);
    yv.push_back (sv[i]);
  }    
  acc = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_cspline, xv.size());
  gsl_spline_init (spline, xv.data(), yv.data(), xv.size());
}
Spline::~Spline ()
{
  if (acc != 0) {
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }
}


int width, height;
size_t snapidx;
vector<array<double, 64> > snaps;
array<double, 64> mean;

int mag = 1;

static cairo_surface_t *surface = NULL;


static void
clear_surface (void)
{
  cairo_t *cr;

  cr = cairo_create (surface);

  cairo_set_source_rgb (cr, 1, 1, 1);
  cairo_paint (cr);

  cairo_destroy (cr);
}

void
get_data (string chanfile, string spikefile, int trigger)
{
  ifstream f(spikefile);
  if (!f)
    error (1, errno, "error opening %s", spikefile.c_str());
  ifstream chan(chanfile);
  if (!chan)
    error (1, errno, "error opening %s", chanfile.c_str());
  string s;
  
  getline (f, s);
  getline (f, s);
  if (s.size() != 13 && s.size() < 15)
    error (1, 0, "bad line length in spike file: %zd", s.size());
  double ticks_per_s = (s.size() == 13 ? 2000 : 10000);
  while (getline (f, s)) {
    int id = atoi (s.substr(0,5).c_str());
    if (id == trigger) {
      double t = atof (s.c_str() + 5);
      double tsamp = t / ticks_per_s * 25000;
      int left = floor (tsamp);
      double delta = tsamp - left;
      if (left < 20) continue;
      array<short, 65> snap65;;
      chan.seekg ((left - 20) * 2);
      chan.read ((char *)snap65.data(), streamsize (snap65.size() * sizeof (short)));
      Spline s (snap65);
      array<double, 64> snap;;
      for (size_t i = 0; i < snap.size(); ++i)
        mean[i] += (snap[i] = s.val (i + delta));
      snaps.push_back (snap);
    }
  }
  f.close ();
  for (auto &m : mean) m /= snaps.size();
  if (snaps.size() == 0) {
    cout << endl;
    error (1, 0, "Code %d not seen in %s\n", trigger, spikefile.c_str());
  }
  
}

static gboolean
calc (gpointer udata)
{
  cairo_t *cr = cairo_create (surface);

  array<double, 64> snap = snaps[snapidx];
  
  cairo_move_to (cr, 0, (-mag * snap[0] + 32768.) / 65536 * height);
  for (size_t i = 1; i < snap.size(); ++i)
    cairo_line_to (cr, 1. * i / snap.size() * width,
                   (-mag * snap[i] + 32768.) / 65536 * height);
  cairo_stroke (cr);
  bool again = ++snapidx < snaps.size();
  if (!again) {
    cairo_move_to (cr, 0, (-mag * mean[0] + 32768.) / 65536 * height);
    for (size_t i = 1; i < mean.size(); ++i)
      cairo_line_to (cr, 1. * i / mean.size() * width,
                     (-mag * mean[i] + 32768.) / 65536 * height);
    cairo_set_source_rgb (cr, 1, 1, 1);
    cairo_stroke (cr);
    
    cairo_move_to (cr, 20. / mean.size() * width, 0);
    cairo_line_to (cr, 20. / mean.size() * width, height);
    cairo_set_source_rgb (cr, 1, 0, 0);
    cairo_stroke (cr);
    cairo_set_source_rgb (cr, 0, 0, 0);

  }
  return again;
}

/* Create a new surface of the appropriate size to store our scribbles */
static gboolean
configure_event_cb (GtkWidget         *widget,
                    GdkEventConfigure *event,
                    gpointer           data)
{
  if (surface)
    cairo_surface_destroy (surface);

  width = gtk_widget_get_allocated_width (widget);
  height = gtk_widget_get_allocated_height (widget);

  surface = gdk_window_create_similar_surface (gtk_widget_get_window (widget),
                                               CAIRO_CONTENT_COLOR,
                                               width, height);

  /* Initialize the surface to white */
  clear_surface ();

  GtkWidget *da = GTK_WIDGET (data);
  if (snapidx >= snaps.size())
    g_idle_add (calc, da);
  snapidx = 0;
  
  /* We've handled the configure event, no need for further processing. */
  return TRUE;
}

/* Redraw the screen from the surface. Note that the ::draw
 * signal receives a ready-to-be-used cairo_t that is already
 * clipped to only draw the exposed areas of the widget
 */
static gboolean
draw_cb (GtkWidget *widget,
         cairo_t   *cr,
         gpointer   udata)
{
  cairo_set_source_surface (cr, surface, 0, 0);
  cairo_paint (cr);

  return FALSE;
}

static gboolean
update (gpointer udata)
{
  GtkWidget *da = GTK_WIDGET (udata);
  gtk_widget_queue_draw (da);
  return true;
}

static gboolean
key_press_event_cb (GtkWidget      *widget,
                    GdkEventKey *event,
                    gpointer        data)
{
  if (event->type != GDK_KEY_PRESS)
    return FALSE;
  guint k = event->keyval;
  if (k == 'q')
    exit (0);
  if (k == '+' || k == '=' || k == '-') {
    if (k == '-') {
      if (mag > 1)
        mag -= 1;
    }
    else mag += 1;
    clear_surface ();
    GtkWidget *da = GTK_WIDGET (data);
    if (snapidx >= snaps.size())
      g_idle_add (calc, da);
    snapidx = 0;
  }  
  return true;
}

static gboolean
button_press_event_cb (GtkWidget      *widget,
                       GdkEventButton *event,
                       gpointer        data)
{
  /* paranoia check, in case we haven't gotten a configure event */
  if (surface == NULL)
    return FALSE;

  if (event->button == GDK_BUTTON_PRIMARY) {
  }
  else if (event->button == GDK_BUTTON_SECONDARY)
    clear_surface ();

  /* We've handled the event, stop processing */
  return TRUE;
}

static void
close_window (void)
{
  if (surface)
    cairo_surface_destroy (surface);

  gtk_main_quit ();
}

static void
activate (void)
{
  GtkWidget *window;
  GtkWidget *frame;
  GtkWidget *drawing_area;

  window = gtk_window_new(GTK_WINDOW_TOPLEVEL);

  gtk_window_set_title (GTK_WINDOW (window), "Drawing Area");

  g_signal_connect (window, "destroy", G_CALLBACK (close_window), NULL);

  gtk_container_set_border_width (GTK_CONTAINER (window), 8);

  frame = gtk_frame_new (NULL);
  gtk_frame_set_shadow_type (GTK_FRAME (frame), GTK_SHADOW_IN);
  gtk_container_add (GTK_CONTAINER (window), frame);

  drawing_area = gtk_drawing_area_new ();
  /* set a minimum size */
  gtk_widget_set_size_request (drawing_area, 300, 300);

  gtk_container_add (GTK_CONTAINER (frame), drawing_area);

  /* Signals used to handle the backing surface */
  g_signal_connect (drawing_area, "draw",
                    G_CALLBACK (draw_cb), NULL);
  g_signal_connect (drawing_area,"configure-event",
                    G_CALLBACK (configure_event_cb), NULL);

  /* Event signals */
  g_signal_connect (drawing_area, "button-press-event",
                    G_CALLBACK (button_press_event_cb), NULL);
  g_signal_connect (window, "key-press-event",
                    G_CALLBACK (key_press_event_cb), NULL);

  /* Ask to receive events the drawing area doesn't normally
   * subscribe to. In particular, we need to ask for the
   * button press and motion notify events that want to handle.
   */
  gtk_widget_set_events (drawing_area, gtk_widget_get_events (drawing_area)
                                     | GDK_BUTTON_PRESS_MASK
                                     | GDK_KEY_PRESS_MASK);

  gtk_widget_show_all (window);

  g_idle_add (calc, drawing_area);
  g_timeout_add (40, update, drawing_area);
}

int
main (int argc, char **argv)
{
  if (argc != 4) {
    cout << "\n  usage: align CHANFILE SPIKEFILE TRIGGER\n\n";
    exit (0);
  }
  get_data (argv[1], argv[2], atoi (argv[3]));

  gtk_init(&argc, &argv);
  
  activate();

  gtk_main();
  
  return 0;
}
