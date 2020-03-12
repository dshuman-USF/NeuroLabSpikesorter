#!/usr/bin/perl

use Spreadsheet::WriteExcel;
($datadir, $prefix, $outdir) = @ARGV;
chomp ($ppath = `which $0`);
chomp ($wd = `pwd`);

if (@ARGV < 2) {
    print "\n";
    print "usage: $0 datadir prefix\n";
    print "\n";
    print "$0 will look for a coordinate spreadsheet file at /oberon/experiments/prefix/prefix.xls\n";
    print "$0 will look for a .nam file at datadir/prefix.nam\n";
    print "$0 will  write  an info file at /oberon/experiments/prefix/prefix_info_orig.xls\n\n";
    exit (0);
}

if (@ARGV > 3) {
    print "\n";
    print "usage: $0 datadir prefix outdir\n";
    print "\n";
    print "$0 will look for a coordinate spreadsheet file at /oberon/experiments/$prefix/$prefix.xls\n";
    print "$0 will look for a .nam file at $datadir/$prefix.nam\n";
    print "$0 will  write  an info file at $outdir/${prefix}_info_orig.xls\n";
    exit (0);
}

$outdir = "/oberon/experiments/$prefix" unless $outdir;

$infocurr = "$outdir/$prefix" . "_info_curr.xls";
if (-e $infocurr) {
    print "\n";
    print "	$infocurr exists, backing it up\n";
    print "	You will need to make sure that any hand-entered edits from the highest-numbered backup\n";
    print "	of $infocurr have been copied\n";
    print "	into the newly created $infocurr.\n\n";
    system "cp -f --backup=numbered $infocurr $infocurr";

    my $parser   = Spreadsheet::ParseExcel->new();
    my $workbook = $parser->parse($infocurr);

    if ( !defined $workbook ) {
        die $parser->error(), ".\n";
    }

    my $worksheet = $workbook->{Worksheet}[0];
    my ( $row_min, $row_max ) = $worksheet->row_range();


    for my $row ( 4 .. $row_max ) {
        my $mchan = $worksheet->get_cell ($row, 1)->value ();
        
        for my $col ( 7 .. 27 ) {
            my $val = $worksheet->get_cell ($row, $col)->value ();
            push @{$aa_sta{$mchan}}, $val;
        }
    }    
}


my $workbook = Spreadsheet::WriteExcel->new("$outdir/$prefix" . "_info_orig.xls");
$worksheet = $workbook->add_worksheet('Info File');
$worksheet2 = $workbook->add_worksheet('AA & STA codes');
$worksheet3 = $workbook->add_worksheet('About Info File');

$worksheet->set_zoom(75);

$mfmt = $workbook->addformat(valign  => 'vcenter', align => 'right', bold => 1, size => 12, bottom => 2);
$worksheet->merge_range('A1:C1', 'EXPERIMENT DATE:', $mfmt);

$worksheet->merge_range('F1:H1', 'EXPERIMENT:', $mfmt);
$worksheet->write('I1', '', $workbook->addformat(valign  => 'vcenter', align => 'center', bold => 0, size => 12,
						 bg_color => 'yellow', bottom => 2, right => 2));

$worksheet->merge_range('J1:L1', 'RECORDING:', $mfmt);
$worksheet->write('M1', '', $workbook->addformat(valign  => 'vcenter', align => 'center', bold => 0, size => 12,
						 bg_color => 'yellow', bottom => 2, right => 2));

$bot2fmt = $workbook->addformat (bottom => 2);
for $col (13..27) { $worksheet->write_blank (0, $col, $bot2fmt); }
$worksheet->write_blank (0, 27, $workbook->addformat (bottom => 2, right => 2));


@hdr_a = ("UNIT NAME","MERGED CH. #","A/P","R/L","DEPTH","dchan","ref");

@hdr_b = ("cord","RLN","vagus","cVRG", "rVRG","rt. PRG","lt. PRG","AA 1 (extra)",
	  "AA 2 (extra)","AA 3 (extra)","phrenic", "RLN","c. vagus","lumbar","cerv. sym.",
	  "exp. laryn","splanchnic","STA 1 (extra)", "STA 2 (extra)","STA 3 (extra)");

@hdr_c = ("COMMENTS (no commas, please)");

$rot = 60;
$hdrfmt_a = $workbook->add_format(bold => 1, size => 12, text_wrap => 1, align => "center");

$hdrfmt_b = $workbook->add_format(bold => 1, size => 12, text_wrap => 0, align => "center", rotation => $rot, left => 1, right => 1, top => 2);
$hdrfmt_b1 = $workbook->add_format(bold => 1, size => 12, text_wrap => 0, align => "center", rotation => $rot, left => 6, top => 2);
$hdrfmt_b2 = $workbook->add_format(bold => 1, size => 12, text_wrap => 0, align => "center", rotation => $rot, right => 6, top => 2);


#$hdrfmt_b = $workbook->add_format(bold => 1, size => 12, text_wrap => 0, align => "center", rotation => $rot, diag_type => 2, top => 2);
#$hdrfmt_b1 = $workbook->add_format(bold => 1, size => 12, text_wrap => 0, align => "center", rotation => $rot, diag_type => 2, top => 2);
#$hdrfmt_b2 = $workbook->add_format(bold => 1, size => 12, text_wrap => 0, align => "center", rotation => $rot, diag_type => 2, top => 2);

$col = 0;

@w = (9.33, 10.22, 9.33, 9.33, 9.33, 9.33, 9.33, 5.89, 5.89, 5.89, 5.89, 5.89,
      5.89, 5.89, 2.78, 2.89, 2.89, 5.89, 5.89, 5.89, 5.89, 5.89, 5.89,
      5.89, 2.78, 2.89, 2.89, 46.67);

$cellfmt = $workbook->add_format(border => 1);


for (@hdr_a) {$worksheet->write(3, $col++, $_, $hdrfmt_a);}
for (@hdr_b) {$worksheet->write(3, $col++, $_, $hdrfmt_b);}
for $n (7,17)  {$worksheet->write(3, $n, $hdr_b[$n-7], $hdrfmt_b1);}
$worksheet->write(3, 26, $hdr_b[26-7], $hdrfmt_b2);
for (@hdr_c) {$worksheet->write(3, $col++, $_, $hdrfmt_a);}

$col = 0;
for (@w) { $worksheet->set_column ($col, $col, $_); $col++ }

$worksheet->set_row (0, 23.4);
$worksheet->set_row (2, 22.8);
$worksheet->set_row (3, 73.2);
for $n (4..200) {$worksheet->set_row ($n, 21.0)}

for $row (4..200) {for $col (0..26) {$worksheet->write_blank ($row, $col, $cellfmt);}}

$cellfmt2 = $workbook->add_format(top => 1, bottom => 1, left => 6);
for $row (4..200) {for $col (7,17,27) {$worksheet->write_blank ($row, $col, $cellfmt2);}}

$aafmt = $workbook->addformat(valign  => 'vcenter', align => 'center', bold => 1, size => 14);
$worksheet->merge_range('J3:Q3', 'AA', $aafmt);

$stafmt = $workbook->addformat(align => 'center', bold => 1, size => 14);
$worksheet->merge_range('T3:AA3', 'STA', $stafmt);

$sta_codes = [
	      "0 -- nothing (no response)",
	      "1 -- peak right - sharp",
	      "       M -- peak right - sharp --> suggest motoneuron",
	      "       P -- peak right - sharp --> suggest premotoneuron",
	      "2 -- peak right - broad",
	      "3 -- peak center",
	      "       D -- peak center, followed by dip",
	      "4 -- trough right",
	      "5 -- trough center",
	      "6 -- multiple peaks and troughs",
	      "7 -- other",
	      "8 -- not tested using this site",
	      "9 -- cardiac EKG"
	      ];

$sz10fmt = $workbook->addformat (size => 10);
$worksheet2->write ('B3', [$sta_codes], $sz10fmt);
$worksheet2->set_column ('A:A', 3.22);
$worksheet2->set_column ('B:B', 51.11);

$aa_codes = [
	     "aa -- cell is antidromically activated from this site",
	     "na -- cell is NOT antidromically activated from this site",
	     "   -- cell was not tested using this site"
	     ];

$worksheet2->write ('B19', [$aa_codes], $sz10fmt);

$bold10fmt = $workbook->addformat (bold => 1, size => 10);
$worksheet2->write ('A1', "STA responses are coded as follows:", $bold10fmt);
$worksheet2->write ('A17', "AA information is coded as follows:", $bold10fmt);

@instr = (
	  "ABOUT INFO FILE:",
	  "",
	  "This file is generated by the program legacy_info_file.pl for legacy experiments and spikesorter for new experiments",
	  "    -legacy_info_file.pl is stored on //Oberon/Experiments/legacy_info_file",
	  "",
	  "The program reads experiment excel files stored in a subdirectory bearing the appropriate experiment name on //Oberon/Experiments",
	  "",
	  "The program generates an excel file called *_info_orig.xls",
	  "",
	  "This file is used to read into xanalysis and the atlas program for future processing",
	  "",
	
	  );

$worksheet3->write ('A3', [["1", "", "", "2", "", "3", "", "4", "",]], $workbook->addformat (align => 'center'));
$worksheet3->write ('B1', [\@instr]);
$worksheet3->set_column ('A:A', 4.11);
$worksheet3->set_column ('B:B', 109.56);


use Spreadsheet::ParseExcel;
use vars '$map';

$infile = "/oberon/experiments/$prefix/$prefix.xls";
if (! -e $infile) {
    $infile = "/oberon/docs/coords/$prefix.xls";
}
$_ = $infile;
if     (/^\.\//) { $ipath = $wd . substr $_, 1; }
elsif (/^\//  )  { $ipath = $_                ; }
else             { $ipath = "$wd/$_"          ; }

$book = Spreadsheet::ParseExcel::Workbook->Parse($infile);

sub check
{
    my ($x, $a,$b,$c,$d,$e,$f,$g,$h) = @_;

    $retval = $sheet->{Cells}[10][ 2+$x]{Val} eq "$b" # DIGITAL CHANNEL
        &&    $sheet->{Cells}[10][ 8+$x]{Val} eq "$c" # A/P
        &&    $sheet->{Cells}[10][ 9+$x]{Val} eq "$d" # R/L
        &&    $sheet->{Cells}[10][10+$x]{Val} eq "$e" # DEPTH
        &&    $sheet->{Cells}[ 2][ 0+$x]{Val} eq "$f" # EXPERIMENT
        &&    $sheet->{Cells}[ 4][ 0+$x]{Val} eq "$g" # RECORDING
        &&    $sheet->{Cells}[ 6][ 0+$x]{Val} eq "$h"; # DATE
    return $retval;
}

sub getvals
{
    $x = shift;
    my $ref_dchan = $sheet->{Cells}[47][3]{Val};
    $refs{$ref_dchan}++;
    $ap_r[$ref_dchan] = $sheet->{Cells}[49][3]{Val};
    $rl_r[$ref_dchan] = $sheet->{Cells}[50][3]{Val};
    $dp_r[$ref_dchan] = $sheet->{Cells}[51][3]{Val};

    for ($r = 11; $sheet->{Cells}[$r][2+$x]{Format}{Fill}[0] == 1; $r++) {
#	print "reading row $r\n";
	$dchan = $sheet->{Cells}[$r][ 2+$x]{Val};
	(exists $ap[$dchan] || exists $rl[$dchan] || exists $dp[$dchan]) && die;
	$ap[$dchan] = $sheet->{Cells}[$r][ 8+$x]{Val};
	$rl[$dchan] = $sheet->{Cells}[$r][ 9+$x]{Val};
	$dp[$dchan] = $sheet->{Cells}[$r][10+$x]{Val};
	$rf[$dchan] = $ref_dchan;
    }
}

$max_ref_mchan = 50;

for (`cat $datadir/$prefix.nam`) {
    chomp;
    ($dchan, $name, $mchan) = split;
    $mchan > $max_ref_mchan || die;
    $dchan{$mchan} = $dchan;
    $name{$mchan} = $name;
}

@vals   = ("DM","DIGITAL CHANNEL","A/P","R/L","DEPTH","EXPERIMENT:","RECORDING:","DATE:");
@novals = (""  ,""               ,""   ,""   ,""     ,""           ,""          ,""     );
@oneval = (""  ,""               ,""   ,""   ,""     ,"(+ = rostral to obex; - = caudal to obex)",
                                                                    ""          ,""     );

sub emsg
{
    $msg = ("The format of\n$ipath\nis not as expected by\n$ppath\n"
	    . "Fix and re-run merge\n");
    system "alert \"$msg\"";
    return 1;
}

for ($n = 0; $n < $book->{SheetCount}; $n++) {
    next if $book->{Worksheet}[$n]{Cells}[0][0]{Val} =~ /^Calibration/;
    next if $book->{Worksheet}[$n]{Cells}[0][0]{Val} =~ /^Using/;
    $sheet = $book->{Worksheet}[$n];
    check (0, @vals) || emsg() && die;
    ($twocol = check (15, @vals)) || check (15, @novals) || check (15, @oneval) || emsg() && die;
    
    getvals (0);
    getvals (15) if $twocol;
}

sub array_letter
{
    my $ref = shift;
    my ($mchan, $dchan);
    keys %dchan;                # reset %dchan iterator
    while (($mchan, $dchan) = each %dchan) {
	if ($rf[$dchan] == $ref) {
	    return substr ($name{$mchan}, 0, 1);
	}
    }
    print "no cell has $ref as a ref\n";
}

$ref_mchan = 0;
for $ref_dchan (sort {$a<=>$b} keys %refs) {
    $ref_mchan++;
    $ref_mchan <= $max_ref_mchan || die;
    $name{$ref_mchan} = array_letter ($ref_dchan) . '0';
    $dchan{$ref_mchan} = $ref_dchan;
}

$worksheet->merge_range('D1:E1', "$sheet->{Cells}[6][2]{Val}",
			$workbook->addformat(valign  => 'vcenter', align => 'center', bold => 0,
					     num_format => 'dd-mmm-yyyy',
					     size => 12, bg_color => 'yellow', bottom => 2, right => 2));

$mfmt2 = $workbook->addformat(valign  => 'vcenter', align => 'center', bold => 0,
			      size => 12, bg_color => 'yellow', bottom => 2, right => 2);

$worksheet->write('I1', "$sheet->{Cells}[2][2]{Val}", $mfmt2);

my $tmp = sprintf "%02d", $sheet->{Cells}[4][2]{Val};
$worksheet->write('M1', "$tmp", $mfmt2);

$numfmt = $workbook->add_format(border => 1, num_format => '0.00');
$numfmt2 = $workbook->add_format(border => 1, num_format => '0');

$orow = 4;
for $mchan (sort {$a<=>$b} keys %dchan) {
    $dchan = $dchan{$mchan};
    $worksheet->write                               ($orow,   0, [$name{$mchan}, $mchan                      ], $cellfmt);
    if (%aa_sta) {
        $worksheet->write ($orow,  7, \@{$aa_sta{$mchan}}, $cellfmt);
    }
    if ($mchan > $max_ref_mchan) {$worksheet->write ($orow,   2, [$ap  [$dchan], $rl  [$dchan], $dp  [$dchan]], $numfmt );}
    else                         {$worksheet->write ($orow,   2, [$ap_r[$dchan], $rl_r[$dchan], $dp_r[$dchan]], $numfmt );}
    $worksheet->write                               ($orow++, 5, [$dchan,        $rf[$dchan]                 ], $numfmt2);
}

$workbook->close();

system ("cp $outdir/${prefix}_info_orig.xls $infocurr")
