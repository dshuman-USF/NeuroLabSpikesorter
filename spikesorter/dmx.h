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


/* dmx.h */


/* added lss Feb5 2008 in an attempt at mac osx porting */
#ifdef __APPLE__
#define fseeko64 fseeko
#define fopen64 fopen
#define ftello64 ftello
typedef long long off64_t;
#endif

typedef struct {
  int type;
  int board;
  int channel;
  int start;
  int end;
  int vscale;
  int voffset;
  int grid_enable;
  int eng_unit_enable;
  int color;
  char unused[32];
} DisplayInfo;
  
  
typedef struct {
  char fid[8];
  int version;
  int board_count;
  int channel_count;
  int swipe_count;
  int display_count;
  DisplayInfo display_info[8];
  short year;
  short month;
  short wday;
  short mday;
  short hour;
  short minute;
  short second;
  short msec;
  int note_size;
  char note[1024];
  char unused[256];
} GlobalHeader;

typedef struct {
  int board;
  int channel;
  int bandwidth;
  int sample_rate;
  int scalar;
  float scale;
  char unit[8];
  char label[32];
  char comment[80];
  int enable;
  float eng_unit_offset;
  int zero_enable;
  float zero_value;
  char unused[40];
} ChannelInfo;  

typedef struct {
  int serial;
  int channel_count;
  ChannelInfo chan[8];
} BoardHeader;
