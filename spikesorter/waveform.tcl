#!/usr/local/bin/waveform

#  Copyright 2005-2020 Kendall F. Morris
#
#   This file is part of the Spiksorter software suite.
#
#   The Spikesorter software suite is free software: you can redistribute
#   it and/or modify it under the terms of the GNU General Public
#   License as published by the Free Software Foundation, either
#   version 3 of the License, or (at your option) any later version.
#
#   The suite is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with the suite.  If not, see <https://www.gnu.org/licenses/>.



if {$argc > 1} {set sample [lindex $argv 1].0} else { set sample 0.0 }
if {$argc > 2} {
    set right_sample [lindex $argv 2].0
} else {
    set right_sample 0
}
set xscale 1.0

set width 1265
set samples_per_width [expr $width / $xscale]
set time_left [format "%.5f" [expr $sample / 25000.0]]
set pick_x -1
set nopaint 0
set ymag 1

proc set_dr_times {} {
    global time_left time_diff time_right samples_per_width sample sample_count width xscale right_sample
    set width [winfo width .sq]
    if {$right_sample > 0 && $right_sample > $sample && $width > 200} {
	set xscale [expr $width / ($right_sample - $sample)]
	set right_sample 0
    }
    set samples_per_width [expr $width / $xscale]
    .sb set [expr $sample / $sample_count]  [expr ($sample + $samples_per_width) / $sample_count] 
    set time_right [format "%.5f" [expr $time_left + $samples_per_width / 25000.0]]
    set time_diff [format "%.5f" [expr $time_right - $time_left]]
    update  idletasks
}

square .sq -type 4 -bg white -xscale $xscale -time $sample -bd 0

frame .sq.fcm -bg red
frame .sq.fym -bg red

frame .fl
label .lleft -textvariable time_left
label .lright -textvariable time_right
label .ldiff -textvariable time_diff
label .lmtim -textvariable mouse_time
label .loff -text " + .00246 = "
label .lmtimo -textvariable mouse_time_offset
label .ly -textvariable cury
label .lmsamp -textvariable mouse_samp
label .lymag -textvariable ymag -fg red

pack .lleft -in .fl -side left
pack .lmtim -in .fl -side left
pack .loff -in .fl -side left
pack .lmtimo -in .fl -side left
pack .ly -in .fl -side left
pack .lmsamp -in .fl -side left
pack .ldiff -in .fl -side left -expand 1
pack .lymag -in .fl -side left
pack .lright -in .fl -side right
pack .fl -fill x

pack .sq -expand 1 -fill both

scrollbar .sb -orient horizontal -command scroll
pack .sb -fill x

if {$argc > 0} {
    set filename [lindex $argv 0]
} else {
    set filename "/cygdrive/d/New_Data_Entry/cygwin/spikesorter/dmx_split_dir/march7_truncated12.chan"
}
wm title . [file tail $filename]

set chan [open_chan $filename]
set sample_count [get_sample_count $chan]

.sq configure -data $chan

proc new_ymag {val} {
    global ymag
    if {$val == 0} {
	set ymag 1
    } elseif {$ymag + $val >= 1} {
	set ymag [expr $ymag + $val];
    }
    .sq configure -yscale [expr .005 * $ymag]
}

bind .sq <Motion> {mouse %W %X %Y}
bind .sq.fcm <Motion> {mouse %W %X %Y}
bind .sq.fym <Motion> {mouse %W %X %Y}

bind .sq <1> {pick %W %X}
bind .sq.fcm <1> {pick %W %X}
bind .sq.fym <1> {pick %W %X}

bind .sq <Shift-1> {unpick %W %X}
bind .sq.fcm <Shift-1> {unpick %W %X}
bind .sq.fym <Shift-1> {unpick %W %X}

bind all <ButtonRelease> {set nopaint 0}
#bind all <KeyPress> {puts %K}
bind all <equal> {puts "$sample $xscale"}
set hold 0
bind all <KeyPress-h> {
    set hold [expr $hold ^ 1]
    if {$hold} {
        .sq configure -bg "sky blue"
    } else {
        .sq configure -bg white
    }
}
bind all <KeyPress-p> {.sq configure -run 1}
bind all <KeyPress-s> {.sq configure -run 0}

bind all <Key-plus>  {new_ymag 1}
bind all <Key-equal>  {new_ymag 1}
bind all <Key-minus> {new_ymag -1}
bind all <Key-0> {new_ymag 0}

set configure_count 0
bind .sq <Configure> {
    if {$configure_count == 0} {
        incr configure_count
        set width 1265
        wm geometry . ${width}x480+0+0
        update idletasks
    } else {
        if {$hold} {
            set xscale [expr %w.0 / $width * $xscale]
        }
        set width %w
    }
    set samples_per_width [expr $width / $xscale]
    set_dr_times
    .sq configure -xscale $xscale
}

bind all <q> exit
bind all <Key-Up> {
    set xscale [expr $xscale / 2.0]
    set samples_per_width [expr $width / $xscale]
    set_dr_times
    .sq configure -xscale $xscale
}
bind all <Key-Down> {
    set xscale [expr $xscale * 2.0]
    set samples_per_width [expr $width / $xscale]
    set_dr_times
    .sq configure -xscale $xscale
}

bind all <Key-Left> {
    set sample [expr $sample - $samples_per_width]
    if {$sample < 0} {set sample 0}
    set time_left [format "%%.5f" [expr $sample / 25000.0]]
    set_dr_times
    .sq configure -time $sample
}
bind all <Key-Right> {
    set sample [expr $sample + $samples_per_width]
    if {$sample < 0 || $sample > $sample_count} {set sample $sample_count }
    set time_left [format "%%.5f" [expr $sample / 25000.0]]
    set_dr_times
    .sq configure -time $sample
}

proc scroll {cmd args} {
    global sample sample_count xscale time_left samples_per_width nopaint
    set amt [lindex $args 0]
    if {$cmd == "moveto"} {
	set sample [expr $amt * $sample_count]
	set nopaint 1
    } elseif {[lindex $args 1] == "units"} {
	set samples_per_pixel [expr 1 / $xscale]
	if {$samples_per_pixel < 1} {set samples_per_pixel 1}
	if {$samples_per_pixel > 1} {set samples_per_pixel [expr $samples_per_pixel * 4]}
	set sample [expr $sample + $amt * $samples_per_pixel]
    } else {
	set sample [expr $sample + $amt * $samples_per_width]
    }
    if {$sample < 0} {set sample 0}
    if {$sample > $sample_count} {set sample $sample_count }
    set time_left [format "%.5f" [expr $sample / 25000.0]]
    set_dr_times
    .sq configure -time $sample
}

frame .sq.fc0 -bg black

proc unpick {cwin cx} {
    global pick_x
    set pick_x -1
    place forget .sq.fc0
}

proc pick {cwin cx} {
    global pick_x sample xscale samples_per_width time_left sample_count
    set cx [expr $cx - [winfo rootx .sq]]
    if {$pick_x == -1} {
	set pick_x $cx
	place .sq.fc0 -in .sq -x $cx -y 0 -width 1 -relheight 1.0
    } else {
	if {$pick_x < $cx} {
	    set x0 $pick_x
	    set x1 $cx
	} else {
	    set x0 $cx
	    set x1 $pick_x
	}		
	set width [winfo width .sq]

	set sampl0 [expr $sample + $x0 / $xscale]
	set sampl1 [expr $sample + $x1 / $xscale]
	set samples_per_width [expr $sampl1 - $sampl0]
	set sample $sampl0
	set xscale [expr $width / $samples_per_width]
	set time_left [format "%.5f" [expr $sample / 25000.0]]

	set_dr_times
	.sq configure -xscale $xscale -time $sample

	set pick_x -1
	place forget .sq.fc0
    }
}

proc mouse {cwin cx cy} {
    global mouse_time mouse_time_offset sample xscale cury mouse_samp ymag
    set height [winfo height .sq]
    set cx [expr $cx - [winfo rootx .sq]]
    set cy [expr $cy - [winfo rooty .sq]]
    set cury [expr ($height / 2 - $cy) * 200 / $ymag]
    set width [winfo width .sq]
    set mouse_time [format "%.5f" [expr ($sample + $cx / $xscale) / 25000.0]]
    set mouse_samp [format "%.0f" [expr ($sample + $cx / $xscale)]]
    set mouse_time_offset [format "%.5f" [expr $mouse_time + .00246]]
    place forget .sq.fcm
    place .sq.fcm -in .sq -x $cx -y 0 -width 1 -relheight 1.0
    place forget .sq.fym
    place .sq.fym -in .sq -x 0 -y $cy -relwidth 1.0 -height 1
}
