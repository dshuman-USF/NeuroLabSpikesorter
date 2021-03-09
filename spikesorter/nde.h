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


/* nde.h */

#define _LARGEFILE64_SOURCE 1
#define _GNU_SOURCE 1

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#if HAVE_UNISTD_H || !HAVE_CONFIG_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <rfftw.h>
#include <errno.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define CLIPPED_SPECIES 2
#undef MIN_W_R
#define PAD 2
#define MPK_R0 ((int)(.0003 * SF + .5))
#define MPK_R ((int)(.0007 * SF + .5))
#define CLIPVAL 27000
#define PRESAMPLES 20
#define WRADIUS 8
#define WRADIUS0 (WRADIUS-3)
#define SF 25000.0
#define SNAPLEN 64
#define NOMINAL_TP_INTERVAL (SF/5)
#define SD 3
#define NOISE_THRESHOLD 11
#define CLUSTER_REMOVAL_THRESHOLD 11
#define CDFSZ 151
#define CDFMULT 10.0
#define MEG *(1<<20)
#define DIE (fprintf (stderr, "fatal error in %s at line %d \n", __FILE__, __LINE__), die(), 0)

#define TMALLOC(buf, n) ((buf = malloc ((n) * sizeof *(buf)))                                                           \
                         && (malloc_debug == 0 || fprintf (stderr, "%d %d %ld\n", __LINE__, 0, (long)buf))) || DIE

#define TCALLOC(buf, n) ((buf = calloc ((n), sizeof *(buf)))                                                            \
                         && (malloc_debug == 0 || fprintf (stderr, "%d %d %ld\n", __LINE__, 0, (long)buf))) || DIE

#define TREALLOC(buf, n) ((malloc_debug == 0 || fprintf (stderr, "%d %ld", __LINE__, (long)buf))        \
                          && (buf = realloc (buf, (n) * sizeof *(buf)))                 \
                          && (malloc_debug == 0 || fprintf (stderr, " %ld\n", (long)buf))) || DIE

#define TMEMCPY(to, from, n) memcpy ((to), (from), (n)*sizeof *(to))
#define TMEMSET(to, from, n) memset ((to), (from), (n)*sizeof *(to))
#define SZMEMCPY(to, from) memcpy ((to), (from), sizeof (to))
#define RND(x) (floor((x) + .5))
#define PWR(n) (fdat[n] * fdat[n] + fdat[N-n] * fdat[N-n])
#undef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#undef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define NEWMIN(a,b) ((b) < (a) ? ((a) = (b)) : (a))
#define NEWMAX(a,b) ((b) > (a) ? ((a) = (b)) : (a))
#define MOD(n,m) ((n) % (m) < 0 ? (n) % (m) + (m) : (n) % (m))
#define PAUSE (printf("Pause - press any key to continue"),getchar())
#define FREE_FB(p) do{free_fb(p);p=0;}while(0)
#define SGN(x) (((x) < 0) ? -1 : (((x) > 0) ? 1 : 0))
typedef unsigned int PIXEL;
extern PIXEL * FrameBuffer;
extern int WIDTH;            
extern float buf[];
extern double mean;
extern double unzapped_sd;
extern double detection_threshold;
extern int starting_sample, bufsamples;
extern char *file_name;

void BresLine(int Ax, int Ay, int Bx, int By, PIXEL Color);

typedef struct {int old_to_new, new_to_old;} Map;

typedef struct
{
  gsl_interp * interp;
  float *y_float;
  /*@observer@*/  short *y_short;
  int size;
  float *c;
} CSpline;

typedef float SnapValArray[SNAPLEN];
typedef double SnapDblArray[SNAPLEN];
typedef struct SnapVals SnapVals;
typedef struct Snap Snap;
typedef struct SnapList SnapList;
typedef struct SnapListList SnapListList;
typedef struct FBufRender FBufRender;
typedef struct SpikeData SpikeData;

struct SpikeData
{
  int partition;
  double sample;
  double distance;
  float scale;
  unsigned int leave_out:1;
  unsigned int nocenter:1;
  unsigned int overlapped:1;
};

extern SpikeData *spikedata;
extern float **spike_waveform;
extern int spike_count, spikedata_alloc;
 
struct SnapVals
{
  float val[SNAPLEN];
};

struct Snap
{
  int count;                    /* the values are the average of this many snaps */
  int rmcnt;
  int copies;                   /* this many snaps converged to these values */
  int refcnt;
  int partition;
  int prevpart;
  int region;
  int peak;
  double radius;
  double shift;
  double sample;
  double narrowness;
  double distance;
  double dr, dw;
  float *val;                   
  float *raw;
  float *buf;
  float *bufptr;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  unsigned char left, right, part;
  unsigned char left0, right0;
  unsigned char endl, endr;
  int noise:2;
  unsigned int close:1;
  unsigned int polarity:1;
  unsigned int halfmatch:1;
  unsigned int recenter:1;
  unsigned int nocenter:1;
  unsigned int skip:1;
  unsigned int free_val:1;
  unsigned int free_raw:1;
  unsigned int done:1;
  unsigned int changed:1;
};

struct SnapList
{
  int count;
  Snap **snap;
};

struct SnapListList
{
  int count;
  SnapList **sl;
};


struct FBufRender {
        int W;
        int H;
        PIXEL *pfb;
};
//void fbrender (SnapList *sl, PIXEL *fbuffer, FBufRender *winp, int imageW, 
void fbrender (SnapList *sl, FBufRender *winp, int imageW, 
                           int imageH, int x0, int y0, int scale, int color);

FBufRender *create_frame_buffer (int W, int H);
void clear_fb(FBufRender *winp);
void resize_fb (FBufRender *winp, int W, int H);
void free_fb(/*@only@*/FBufRender *winp);

struct segment {
  double loc0;
  unsigned long tp0;
  double m;
};

struct tp_table {
  double *loc;
  unsigned long tp_count;
  unsigned long tp_alloc;
  struct segment *sg;
  unsigned long sg_count;
  short *buf;
  unsigned long s_count;
  unsigned long first_sample;
};

struct bp_table {
  double *loc;
  unsigned long count;
};

#define save_var(name, len)                             \
do {                                                    \
  FILE *f;                                              \
                                                        \
  (f = fopen (#name, "wb")) || DIE;                     \
  fwrite (name, sizeof *name, len, f) == len || DIE;    \
  fclose (f);                                           \
}  while (0)

#define restore_var(name, len)                          \
do {                                                    \
  FILE *f;                                              \
                                                        \
  (f = fopen (#name, "rb")) || DIE;                     \
  filesize (f) == len * sizeof *name || DIE;            \
  fread (name, sizeof *name, len, f) == len || DIE;     \
  fclose (f);                                           \
}  while (0)

#define save_snaps(sl)                          \
do {                                            \
  SnapValArray *sl##_snaps;                     \
                                                \
  sl##_snaps = save_snaplist (sl, 0);           \
  save_var (sl##_snaps, sl->count);             \
  free (sl##_snaps);                            \
} while (0)

#define restore_snaps(sl)                       \
do {                                            \
  SnapValArray *sl##_snaps;                     \
                                                \
  TMALLOC (sl##_snaps, sl->count);              \
  restore_var (sl##_snaps, sl->count);          \
  restore_snaplist (sl, sl##_snaps);            \
  free (sl##_snaps);                            \
} while (0)

#define FREE2D(p, count)                                        \
{                                                               \
  int free2d_idx;                                               \
  for (free2d_idx = 0; free2d_idx < (count); free2d_idx++)      \
    free ((p)[free2d_idx]);                                     \
  free (p);                                                     \
  (p) = 0;                                                      \
}

#define sample_to_sample(tpfrom, tpto, loc)             \
tp_to_sample (tpto, sample_to_tp (tpfrom, loc))

struct tp_table *find_tp_locs (char *filename);
double sample_to_tp (struct tp_table *tp, double loc);
double tp_to_sample (struct tp_table *tp, double tpt);
double corr_filter (fftw_real *tdat, int N, double freq, double tol);
double *median_array (double *val, int count, int ary_sz, double *m);
double *mean_array (double *val, int count, int ary_sz, double *m);
double radius_max (double *val, int size, int radius, int index);
void display_image (FBufRender *p, int x, int y, char *win_name);
float qselect(unsigned long k, unsigned long n, float arr[]);
void ksone(float data[], unsigned long n, float (*func)(float), float *d, float *prob);
double fmin2 (double (*func)(double x, void *param, int debug), void *p, double startx, double radius, double min, double max, int *status);
double fmin_y (void);
void fblines (FBufRender *winp, int imageW, int imageH, int x0, int y0, int color);


/*@dependent@*/FILE *create_pdf_file (char *f);
void pdflines (FILE *f, int imageW, int imageH, int x0, int y0, unsigned color);
void pdfrender (SnapList *sl, FILE *f, int imageW, int imageH, int x0, int y0, int scale, unsigned color);
void pdfdisplay (FILE *f);
void start_pdf_page (FILE *f);
void end_pdf_page (FILE *f);
void pdftext (FILE *f, int x, int y, char *txt);
void setcolor (FILE *f, unsigned color);

size_t dmx_read (void *buf, size_t size, size_t count, void *p);
void * dmx_open (char *file_name, size_t bufsize);

#ifdef __CYGWIN__
typedef long long off64_t;
int fseeko64 (void *stream, off64_t offset, int whence);
off64_t ftello64 (void *stream);
int fclose64 (void *stream);
void *fopen64 (const char *filename, const char *mode);
size_t fread64 (void * buf, size_t size, size_t count, void * fp);
#else
#define fread64 fread
#define fclose64 fclose
#endif

void create_ach (SnapList *sl, FILE *f, int x0, int y0, float bin_width, int axis_length, unsigned color);
void create_timeline (SnapList *sl, FILE *f, int length, int x0, int y0);
void create_variance_histograms (SnapList *sl, SnapList *noise, FILE *f, int x0, int y0);
void die (void);

double residue (float *);
float *init_residue (SnapList *sl, int first_snap, int last_snap);
double local_cubic (double x, float *y, int count, double edge);

void draw_snap (float *v, /*@observer@*/ char *color);
void clear (void);
int show (void);
int show_nowait (void);
void window (void);
void wclose (void);
void draw_peak (float *v, /*@observer@*/ char *color);
void draw_text (/*@unique@*/ /*@observer@*/ char *t) ;
void whiten_samples (float *v, int count);
int xor_snap (float *v, int wait) ;
void get_click (double *xp, double *yp);
void add_snap (SnapList *sl, Snap *s);
void whiten_snap (float *v);
SnapList * center_points_full (SnapList *sl, Snap *s);
void free_snaplist (SnapList **slp);
float *optavg (SnapList *sl, int first_snap, int last_snap, float *p0, float *buf, int starting_sample, int bufsamples);
double map_radius (double rin, int len);
void show_centers (SnapList *avg);
char * change_filetype (char *file_name, char *suffix);
double white_distance_from_zero (float *a, int count);
double lp_filter (float *v);
void lpfilter_snap (float *v, float *out) ;
extern double *lp_coeff;
double lp_filter_dbl (double *v) ;
extern double cov_rep_root[SNAPLEN][SNAPLEN];
CSpline * cspline_new (/*@only@*/ float *y_float, /*@observer@*/ short *y_short, int size) ;
double distance (float *a, float *b);
void mem_use (char *file, int line);
double norm_distance_from_zero (float *a);
void note (char *fmt, ...);
void shift_snap (float *v, double shift_amt);

struct malloc_chunk {

  size_t      prev_size;  /* Size of previous chunk (if free).  */
  size_t      size;       /* Size in bytes, including overhead. */

  struct malloc_chunk* fd;         /* double links -- used only if free. */
  struct malloc_chunk* bk;
};

typedef struct malloc_chunk* mfastbinptr;
typedef struct malloc_chunk* mchunkptr;

struct malloc_state {
  size_t  max_fast;
  mfastbinptr      fastbins[10];
  mchunkptr        top;
  mchunkptr        last_remainder;
  mchunkptr        bins[96 * 2];
  unsigned int     binmap[3+1];
  unsigned long     trim_threshold;
  size_t  top_pad;
  size_t  mmap_threshold;
  int              n_mmaps;
  int              n_mmaps_max;
  int              max_n_mmaps;
  unsigned int     pagesize;
  unsigned int     morecore_properties;
  size_t  mmapped_mem;
  size_t  sbrked_mem;
  size_t  max_sbrked_mem;
  size_t  max_mmapped_mem;
  size_t  max_total_mem;
};
extern struct malloc_state av_;

#define SIZE_SZ 4
#define mem2chunk(mem) ((mchunkptr)((char*)(mem) - 2*SIZE_SZ))

#define IS_MMAPPED 0x2
#define chunk_is_mmapped(p) ((p)->size & IS_MMAPPED)

#define PREV_INUSE 0x1
#define SIZE_BITS (PREV_INUSE|IS_MMAPPED)
#define chunksize(p)         ((p)->size & ~(SIZE_BITS))
#define chunk_at_offset(p, s)  ((mchunkptr)(((char*)(p)) + (s)))

int malloc_check (void);
extern int malloc_debug;
void myfree (int line, /*@only@*/ void *p);
extern time_t start_time;
extern bool use_best_attempt_flag;

#if !HAVE_ASPRINTF && HAVE_CONFIG_H
int asprintf(char **buffer, char *fmt, ...);
#endif

double ss_cspline_eval (const gsl_spline *spline, double x, gsl_interp_accel *a);
void ss_spline_free (gsl_spline * spline);
int find_detection_threshold (char *file_name) ;
