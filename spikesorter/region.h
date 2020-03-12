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

/* region.h */


#include <gsl/gsl_spline.h>

#define BEFORE (SNAPLEN-1)
#define AFTER (SNAPLEN-1)
#define DATALEN (BEFORE+AFTER+1)
  
typedef struct {
  int sample, count;
  char type;
  unsigned int clipped:1;
} Region;

typedef struct {
  int loc;
  char lfar, lnear, rnear, rfar;
} Peak;
  
typedef struct {
  gsl_interp_accel *acc;
  gsl_spline *spline;
  int pn;
  int spike_count;
  int pkcnt;
  int *clip;
  int clipcnt;
  double *offset;
  /*@dependent@*/ Peak *peak;
  short cluster;
  unsigned int mark:1;
} RegionM;

typedef struct {
  short *data;
  float *fdata;
  Peak *peak;
  Region *region;
  RegionM *regionm;
  int region_count, data_count, peak_count;
  int region_alloc, data_alloc, peak_alloc;
  int *peak_region;
  char *done;
} RegionData;

typedef struct {
  char type;
  /*@dependent@*/ Peak *peak;
  /*@dependent@*/float *cdata;
  float *val;
  /*@dependent@*/float *raw;
  /*@observer@*/ RegionM *regionm;
  int window_size;
  int ccn;
  double scale;
} Seed;

typedef struct {
  gsl_interp_accel *acc;
  gsl_spline *spline;
  Peak peak[3];
  float raw[SNAPLEN];
  int type;
  int peak_count;
} Center;

typedef struct {
  int ccn;
  int peak_count;
  double start;
  double distance;
  double scale;
  float raw[SNAPLEN];
  RegionData *p;
  unsigned int overlapped:1;
} Component;

typedef struct {
  int component_count;
  Component *component;
  double residue;
} Decomposition;

typedef struct {
  int max_partition;
  /*@dependent@*/RegionData *p0;
  /*@owned@*/RegionData *p;
  Center *center;
  int center_count;
  int pn;
  int rn;
  short *region_data_ptr;
  float *region_data_copy;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  int region_data_alloc;
  Seed seed;
  Component *component;
  int component_count;
  int component_alloc;
  double min_distance;
  Decomposition *decomposition;
  int decomposition_count;
  int decomposition_alloc;
  int *clip;
  int clipcnt;
  Map *map;
  int *pkloc;
  int *cpkloc;
  double good_lo, good_hi;
  unsigned int flat:2;
  unsigned int flag:1;
} DecomposeData;

int region_add_data (RegionData *p, int b, int lfar, int lnear, int rnear, int rfar, /*@only@*/ float *buf, int bufsamples, int starting_sample, int preexisting);
int region_check (RegionData *p);
void region_write (RegionData *p, char *filename);
RegionData * region_read (char *filename);
void region_spline_init (RegionData *p);
SnapList * region_snaplist (RegionData *p);
double region_spline_peak (RegionData *p, Snap *s);
void region_clipped (RegionData *p);
void region_show_peaks (RegionData *p, double);
void region_classify (RegionData *p);
SnapList *region_cluster (DecomposeData *dd);
void region_optavg (RegionData *p);
void region_show (RegionData *p, int rns);
void decompose_unclassified_regions (DecomposeData *dd);
void region_info (RegionData *p, int rn);
void region_write_centers (DecomposeData *dd, char *file_name);
SnapList *region_read_centers (DecomposeData *dd, char *file_name);
void region_unclassified (RegionData *p, int);
int revisit_assignments (DecomposeData *dd, SnapList *cc, SnapList *noise);
double one_region_noise (float *p, int count, int *maxidx) ;
int find_mpk_spikes (RegionData *rd, /*@only@*/ float *buf, int bufsamples, int starting_sample, int buf_start, int preexisting);
RegionData * region_read_justone (char *filename, int rn);
void read_spikedata (char *file_name, char *ext);
void region_ctr_init (Center *ctr);
void write_spikedata (char *file_name, int);
