#!/usr/bin/perl
# CONEVA - Contact assessment toolkit
# Updated by Badri Adhikari on Jan 2020 (only calculates precision for one rr file)
# Written by Badri Adhikari and Jackson Nowotny - May 2016

# References:
# http://predictioncenter.org/casproll/index.cgi?page=format#RR (RR file format)
# http://pdg.cnb.uam.es/eva/con/predictionFormats.html (EVAcon format)
# Marks, D.S., et al., Protein 3D structure computed from evolutionary sequence variation. PloS one, 2011. 6(12): p. e28766. (for calculating spread and FP error)
# Monastyrskyy, B., et al., Evaluation of residue–residue contact predictions in CASP9. Proteins: Structure, Function, and Bioinformatics, 2011. 79(S10): p. 119-125. (precision, coverage, Xd)
# Monastyrskyy, B., et al., New encouraging developments in contact prediction: Assessment of the CASP11 results. Proteins, 2015. (for MCC and Jaccard)

use strict;
use warnings;
use Carp;
use File::Basename;
use File::Temp qw(tempfile);
use Getopt::Long;

# User inputs
my ($in_rr, $in_pdb_file, $atom_type, $d_threshold, $min_seq_sep, $max_seq_sep, $jaccard_selct, $flag_show_all);

GetOptions(
	"rr=s"	=> \$in_rr,
	"pdb=s"	=> \$in_pdb_file,
	"atom=s"=> \$atom_type,
	"d=s"	=> \$d_threshold,
	"smin=i"=> \$min_seq_sep,
	"smax=i"=> \$max_seq_sep,
	"all"	=> \$flag_show_all,
	"jacc=s"=> \$jaccard_selct)
or confess "ERROR! Error in command line arguments!";

print_usage() if not defined $in_pdb_file;

# Defaults
$atom_type    = "cb"  if !$atom_type;
$d_threshold  = 8.0   if !$d_threshold;
$max_seq_sep  = 10000 if !$max_seq_sep;
$min_seq_sep  = 24    if !$min_seq_sep;

our %AA3TO1 = qw(ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V UNK -);
our %AA1TO3 = reverse %AA3TO1;

# Precompute PDB stats
my $seq_pdb = seq_chain_with_gaps($in_pdb_file);
my %pdb_rnum_rname = res_num_res_name($in_pdb_file);
my %pdb_dist;
my %pdb_dist_atoms;
my %xyz = xyz_pdb($in_pdb_file);
foreach my $row1(keys %xyz){
	my @row1 = split /\s+/, $row1;
	my $res1 = $row1[0];
	my $atm1 = $row1[1];
	if ($atom_type eq "heavyatoms"){
		next if not ($atm1 eq "N" or $atm1 eq "CA" or $atm1 eq "C" or $atm1 eq "O");
	}
	if ($atom_type eq "ca"){
		next if not $atm1 eq "CA";
	}
	if ($atom_type eq "cb"){
		next if not $atm1 eq return_cb_or_ca_atom($seq_pdb, $res1);
	}
	foreach my $row2(keys %xyz){
		my @row2 = split /\s+/, $row2;
		my $res2 = $row2[0];
		my $atm2 = $row2[1];
		if ($atom_type eq "heavyatoms"){
			next if not ($atm2 eq "N" or $atm2 eq "CA" or $atm2 eq "C" or $atm2 eq "O");
		}
		if ($atom_type eq "ca"){
			next if not $atm2 eq "CA";
		}
		if ($atom_type eq "cb"){
			next if not $atm2 eq return_cb_or_ca_atom($seq_pdb, $res2);
		}
		next if $res1 >= $res2;
		next if abs($res1 - $res2) < $min_seq_sep;
		next if abs($res1 - $res2) > $max_seq_sep;
		my $d = calc_dist($xyz{$row1}, $xyz{$row2});
		if (not defined $pdb_dist{$res1." ".$res2}){
			$pdb_dist{$res1." ".$res2} = $d;
			$pdb_dist_atoms{$res1." ".$res2} = "$atm1 $atm2";
		}
		if ($pdb_dist{$res1." ".$res2} > $d){
			$pdb_dist{$res1." ".$res2} = $d;
			$pdb_dist_atoms{$res1." ".$res2} = "$atm1 $atm2";
		}
	}
}
my %pdb_cont_dist;
foreach (sort keys %pdb_dist){
	$pdb_cont_dist{$_} = $pdb_dist{$_} if $pdb_dist{$_} < $d_threshold;
}

# Load top contact subsets
my %top_count       = ();
my %top_count_order = ();
my $ungapped_seq = seq_chain($in_pdb_file); # use Native's Length as reference
#$top_count{"5"}  = 5;
#$top_count{"L/10"} = int(0.1 * length($ungapped_seq) + 0.5);
$top_count{"L/5"}  = int(0.2 * length($ungapped_seq) + 0.5);
#$top_count{"L/2"}  = int(0.5 * length($ungapped_seq) + 0.5);
$top_count{"L"}    = length($ungapped_seq);
#$top_count{"2L"}   = int(2.0 * length($ungapped_seq));
$top_count{"Nc"}  = scalar keys %pdb_cont_dist;
#$top_count_order{"5"}    = 1;
#$top_count_order{"L/10"} = 2;
$top_count_order{"L/5"}  = 3;
#$top_count_order{"L/2"}  = 4;
$top_count_order{"L"}    = 5;
#$top_count_order{"2L"}   = 6;
$top_count_order{"Nc"}  = 7;

if($jaccard_selct){
	confess "Cannot calculate Jaccard Similarity for $jaccard_selct! Options are 5, L/10, L/5, L/2, L, 2L\n" if not defined $top_count{$jaccard_selct};
}

# Load all contacts to be assessed (including pdb contacts)
my @rr_list = ();
if ($in_rr){
	if (-d $in_rr){
		@rr_list = load_rr($in_rr);
	}
	elsif (-f $in_rr){
		push @rr_list, $in_rr;
	}
}

# Load all rr sequences
my %seq_rr;
foreach my $rr (@rr_list){
	if (rr_has_seq($rr)){
		$seq_rr{$rr} = seq_rr($rr);
	}
	else{
		$seq_rr{$rr} = $seq_pdb;
		$seq_rr{$rr} =~ s/[A-Z]/-/g;
	}
}

# Throw warnings when sequences mismatch
foreach my $rr (@rr_list){
	if ($seq_rr{$rr} ne $seq_pdb){
		warn "\nSequence mismatch!\n";
		warn "".$seq_pdb."[".basename($in_pdb_file)."]\n";
		warn "".$seq_rr{$rr}."[".basename($rr)."]\n";
	}
}

# Load contacts into memory
my %all_rr;
foreach my $rr (@rr_list){
	my %rr_needed = rr_rank_rows_hash($rr, 100000);
	$all_rr{"".basename($rr)." ALL"} = \%rr_needed;
}
foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
	foreach my $rr (@rr_list){
		my %rr_needed = rr_rank_rows_hash($rr, $top_count{$top});
		$all_rr{"".basename($rr)." ".$top} = \%rr_needed;
	}
}
my %all_rr_above_thres_conf;
foreach my $rr (@rr_list){
	my %rr_selected = rr_rank_rows_hash($rr, 100000, 0.5);
	$all_rr_above_thres_conf{"".basename($rr)} = \%rr_selected;
}

# Add native contacts as one of the contacts (for comparison)
my ($pdb_rr_fh, $pdb_rr_file) = tempfile();
open TEMP, ">$pdb_rr_file" or confess $!;
print TEMP "$seq_pdb\n";
foreach (sort keys %pdb_cont_dist){
	print TEMP $_." 0 $d_threshold ".$pdb_cont_dist{$_}."\n";
}
close TEMP;
push @rr_list, basename($in_pdb_file);
foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
	my %rr_needed = rr_rank_rows_hash($pdb_rr_file, $top_count{$top});
	$all_rr{"".basename($in_pdb_file)." ".$top} = \%rr_needed;
}
my %rr_all_true = rr_rank_rows_hash($pdb_rr_file, 100000);
$all_rr_above_thres_conf{"".basename($in_pdb_file)} = \%rr_all_true;
$all_rr{"".basename($in_pdb_file)." ALL"} = \%rr_all_true;
system_cmd("rm -f $pdb_rr_file");

# Contact counts
print "PDB     : $in_pdb_file\n";
print "RR      : $in_rr\n" if $in_rr;
print "L       : ".length(seq_chain($in_pdb_file))." (pdb's chain length without gaps)\n";
print "Nc      : ".(scalar keys %pdb_cont_dist)." (pdb's contact count)\n";
print "Seq Sep : $min_seq_sep to ".(($max_seq_sep == 10000) ? "INF": $max_seq_sep)."\n";
print "PDB-Seq : ".seq_chain_with_gaps($in_pdb_file)."\n";
print "\n";
printf "%-50s", "CONTACT-COUNTS";
foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
	$top = "Top-$top";
	printf "%-10s", $top;
}
print "\n";

foreach my $rr (sort @rr_list){
	printf "%-50s", "".basename($rr)." (count)";
	foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
		my %this_rr = %{$all_rr{"".basename($rr)." ".$top}};
		printf "%-10s", (scalar keys %this_rr);
	}
	print "\n";
}

if ((scalar keys %pdb_cont_dist) < 1){
	print "\nERROR!\n";
	print "Number of true contacts in the pdb file is 0! Quitting..\n";
	exit 1;
}

# Precision
print "\n";
printf "%-50s", "PRECISION";
printf "%-10s", "Top-L/5";
printf "%-10s", "Top-L";
printf "%-10s", "Top-Nc";

print "\n";
foreach my $rr (sort @rr_list){
	printf "%-50s", "".basename($rr)." (precision)";
	my %this_rr = %{$all_rr{"".basename($rr)." L/5"}};
	printf "%-10s", calc_precision(\%this_rr);
	%this_rr = %{$all_rr{"".basename($rr)." L"}};
	printf "%-10s", calc_precision(\%this_rr);
	%this_rr = %{$all_rr{"".basename($rr)." Nc"}};
	printf "%-10s", calc_precision(\%this_rr);
	print "  $in_pdb_file\n";
	last;
}

sub rrhash2dist{
	# Empty input returns empty output
	my $rrhash = shift;
	my %rr_hash = %{$rrhash};
	my %output = ();
	foreach (keys %rr_hash){
		my @C = split /\s+/, $rr_hash{$_};
		next if not defined $pdb_dist{$C[0]." ".$C[1]};
		$output{$rr_hash{$_}} = $pdb_dist{$C[0]." ".$C[1]};
	}
	return %output;
}

sub calc_precision{
	my $rrhash = shift;
	my %rr = %{$rrhash};
	confess "Cannot calculate precision! Empty selected contacts" if not scalar keys %rr;
	my %distances = rrhash2dist(\%rr, \%pdb_dist);
	confess "Distances could not be calculated for selected contacts! Something went wrong!" if not scalar keys %distances;
	my $satisfied = 0;
	foreach (sort {$distances{$a} <=> $distances{$b}}keys %distances){
		my @R = split /\s+/, $_;
		$satisfied++ if $distances{$_} <= $R[3];
	}
	return sprintf "%.2f", 100 * ($satisfied/(scalar keys %distances));
}

sub calc_coverage{
	my $rrhash = shift;
	my %rr = %{$rrhash};
	confess "Cannot calculate precision! Empty selected contacts" if not scalar keys %rr;
	my $total = scalar keys %pdb_cont_dist;
	my $covered = 0;
	my %two_col_rr = ();
	foreach (keys %rr){
		my @R = split /\s+/, $rr{$_};
		$two_col_rr{$R[0]." ".$R[1]} = 1;
	}
	foreach (sort keys %pdb_cont_dist){
		$covered++ if defined $two_col_rr{$_};
	}
	return sprintf "%.2f", 100 * ($covered/$total);
}

sub calc_fp_error{
	# In reference to the EVFOLD paper
	my $rrhash = shift;
	my %rr = %{$rrhash};
	confess "Cannot calculate fp! Empty selected contacts" if not scalar keys %rr;
	my %distances = rrhash2dist(\%rr, \%pdb_dist);
	confess "Distances could not be calculated for selected contacts! Something went wrong!" if not scalar keys %distances;
	my $deviation = 0;
	foreach (keys %rr){
		my @P = split /\s+/, $rr{$_};
		confess "One of the residues in ".$P[0]." ".$P[1]." is not defined in native_pdb!" if not defined $pdb_dist{$P[0]." ".$P[1]};
		my $d = $pdb_dist{$P[0]." ".$P[1]};
		if ($d < $P[3]){
			$deviation += 0.0;
		}
		else{
			$deviation += ($d - $P[3]);
		}
	}
	return sprintf "%.2f", $deviation/(scalar keys %rr);
}

sub calc_xd{
	my $rrhash = shift;
	my $xd_calc = 0;
	my %rr = %{$rrhash};
	my %pred_cont_dist = rrhash2dist(\%rr, \%pdb_dist);
	confess "ERROR! Empty input rr!" if not scalar keys %pred_cont_dist;
	for(my $i = 1; $i <= 15; $i++){
		my $dmin = 4 * ($i - 1);
		my $dmax = 4 * $i;
		my $ppi = count_in_range_xd($dmin, $dmax, \%pred_cont_dist)/(scalar keys %pred_cont_dist);
		my $pai = fraction_of_min_true_dists($dmin, $dmax);
		my $sum = 100 * ($ppi - $pai)/($i);
		$xd_calc += $sum;
	}
	return sprintf "%.2f", $xd_calc;
}

sub count_in_range_xd{
	my $dmin = shift;
	my $dmax = shift;
	my $ref_hash = shift;
	my %cont_dist = %{$ref_hash};
	my $count = 0;
	foreach my $c (keys %cont_dist){
		$count++ if ($cont_dist{$c} >= $dmin and $cont_dist{$c} <= $dmax);
	}
	return $count;
}

sub fraction_of_min_true_dists{
	my $dmin = shift;
	my $dmax = shift;
	my $count = 0;
	foreach (keys %pdb_dist){
		# Keeping < $dmax to match CASP numbers
		$count++ if ($pdb_dist{$_} >= $dmin and $pdb_dist{$_} < $dmax);
	}
	return $count/(scalar keys %pdb_dist);
}


sub calc_spread{
	# Spread: distribution (spread) of the contacts along the chain and
	# over the structure of the protein, by measuring the mean of the distance
	# from every experimental (crystal structure) contact to the nearest predicted contact
	my $rrhash = shift;
	my %rr = %{$rrhash};
	my %two_col_rr = ();
	foreach (keys %rr){
		my @R = split /\s+/, $rr{$_};
		$two_col_rr{$R[0]." ".$R[1]} = 1;
	}
	if (not scalar keys %two_col_rr){
		return "-";
	}
	confess "ERROR! Empty input rr!" if not scalar keys %two_col_rr;
	my $spread = 0;
	confess "There are no contacts from the native!" if not scalar keys %pdb_cont_dist;
	foreach my $nat(keys %pdb_cont_dist){
		my $min = 10000;
		foreach my $pred(keys %two_col_rr){
			my @N = split /\s+/, $nat;
			my @P = split /\s+/, $pred;
			my $d = sqrt(($N[0]-$P[0])*($N[0]-$P[0])+($N[1]-$P[1])*($N[1]-$P[1]));
			$min = $d if $d < $min;
		}
		$spread += $min;
	}
	return sprintf "%.2f", $spread/(scalar keys %pdb_cont_dist);
}

sub calc_mcc{
	my $rrhash = shift;
	my %rr = %{$rrhash};
	if (not scalar keys %rr){
		return "-";
	}
	my %distances = rrhash2dist(\%rr, \%pdb_dist);
	confess "Distances could not be calculated for selected contacts! Something went wrong!" if not scalar keys %distances;
	my $satisfied = 0;
	foreach (sort {$distances{$a} <=> $distances{$b}}keys %distances){
		my @R = split /\s+/, $_;
		$satisfied++ if $distances{$_} <= $R[3];
	}
	my $tp = $satisfied;
	my $fp = (scalar keys %rr) - $satisfied;
	my $tn = (scalar keys %pdb_dist) - (scalar keys %pdb_cont_dist) - $fp;
	my $fn = (scalar keys %pdb_cont_dist) - $satisfied;
	my $mcc = ($tp*$tn - $fp*$fn)/sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn));
	return sprintf "%.2f", $mcc;
}

sub calc_jaccard_similarity{
	my $rr_set1  = shift;
	my $rr_set2  = shift;
	my $NEIGHBOR = shift;
	my %rr1input = %{$rr_set1};
	my %rr2input = %{$rr_set2};
	my %rr1 = ();
	my %rr2 = ();
	my %union = ();
	my %intersect = ();
	foreach (keys %rr1input){
		my @C = split /\s+/, $rr1input{$_};
		$rr1{$C[0]." ".$C[1]} = 1;
	}
	foreach (keys %rr2input){
		my @C = split /\s+/, $rr2input{$_};
		$rr2{$C[0]." ".$C[1]} = 1;
	}
	# add to union set
	foreach (keys %rr1){
		my @C = split /\s+/, $_;
		my $flag_exists = 0;
		for (my $i = -($NEIGHBOR); $i <= $NEIGHBOR; $i++){
			for (my $j = -($NEIGHBOR); $j <= $NEIGHBOR; $j++){
				if (defined $union{($C[0] + $i)." ".($C[1] + $j)}){
					$flag_exists = 1;
					last;
				}
			}
		}
		$union{$_} = 1 if ! $flag_exists;
	}
	foreach (keys %rr2){
		my @C = split /\s+/, $_;
		my $flag_exists = 0;
		for (my $i = -($NEIGHBOR); $i <= $NEIGHBOR; $i++){
			for (my $j = -($NEIGHBOR); $j <= $NEIGHBOR; $j++){
				if (defined $union{($C[0] + $i)." ".($C[1] + $j)}){
					$flag_exists = 1;
					last;
				}
			}
		}
		$union{$_} = 1 if ! $flag_exists;
	}
	# add to intersection set
	foreach (keys %union){
		my @C = split /\s+/, $_;
		my $flag_exists_rr1 = 0;
		my $flag_exists_rr2 = 0;
		for (my $i = -($NEIGHBOR); $i <= $NEIGHBOR; $i++){
			for (my $j = -($NEIGHBOR); $j <= $NEIGHBOR; $j++){
				$flag_exists_rr1 = 1 if defined $rr1{($C[0] + $i)." ".($C[1] + $j)};
				$flag_exists_rr2 = 1 if defined $rr2{($C[0] + $i)." ".($C[1] + $j)};
				last if ($flag_exists_rr1 and $flag_exists_rr2);
			}
		}
		$intersect{$_} = 1 if ($flag_exists_rr1 and $flag_exists_rr2);
	}
	return sprintf "%.2f", (scalar keys %intersect)/((scalar keys %union));
}

sub load_rr{
	my $dir_rr = shift;
	confess ":( directory $dir_rr does not exist!" if not -d $dir_rr;
	my @rr_list = ();
	my @all_files = <$dir_rr/*>;
	foreach my $file (@all_files){
		next if length($file) < 2;
		open FILE, $file or confess $!;
		while(<FILE>){
			chomp $_;
			$_ =~ s/\r//g;
			$_ =~ s/^\s+//;
			next unless $_ =~ /^[0-9]/;
			my @C = split /\s+/, $_;
			confess "column 1 not defined in line [$_] in $file" if not defined $C[0];
			confess "column 2 not defined in line [$_] in $file" if not defined $C[1];
			confess "column 3 not defined in line [$_] in $file" if not defined $C[2];
			confess "column 4 not defined in line [$_] in $file" if not defined $C[3];
			confess "column 5 not defined in line [$_] in $file" if not defined $C[4];
			last;
		}
		close FILE;
		push @rr_list, $file;
	}
	confess "ERROR! Directory $dir_rr has no pdb files!\n" unless(@rr_list);
	return @rr_list;
}

sub rr_rank_rows_hash{
	my $file_rr    = shift;
	my $count      = shift;
	my $conf_thres = shift;
	confess "Input file not defined" if not defined $file_rr;
	confess "Input file $file_rr does not exist!" if not -f $file_rr;
	confess "No contact count!" if not defined $count;
	my ($rr_temp_fh, $rr_temp_file) = tempfile();
	system_cmd("cp $file_rr $rr_temp_file");
	sort_rr_file_by_confidence($rr_temp_file);
	my %rr = ();
	my $i = 1;
	open RR, $rr_temp_file or confess $!;
	while(<RR>){
		my $row = $_;
		next unless $row =~ /[0-9]/;
		chomp $row;
		$row =~ s/\r//g;
		$row =~ s/^\s+//;
		next unless $row =~ /^[0-9]/;
		my @C = split /\s+/, $row ;
		confess "Expecting a pair in row [".$row."]!\n" if (not defined $C[0] || not defined $C[1]);
		confess "Confidence column not defined in row [".$row."] in file $file_rr!\nPlease make sure that the input RR file is in 5-column format!" if not defined $C[4];
		# Skip if native does not have the residue
		next if not defined $pdb_rnum_rname{$C[0]};
		next if not defined $pdb_rnum_rname{$C[1]};
		# Smaller residue number should come first
		$row = $C[1]." ".$C[0]." ".$C[2]." ".$C[3]." ".$C[4] if $C[0] > $C[1];
		if (defined $conf_thres){
			next if $C[4] < $conf_thres;
		}
		# Select only LR, MR, SR or all
		my $d = abs($C[0]-$C[1]);
		next if $d < $min_seq_sep;
		next if $d > $max_seq_sep;
		$rr{$i} = $row;
		$i++;
		last if $i > $count;
	}
	close RR;
	system_cmd("rm -f $rr_temp_file");
	return %rr;
}

sub sort_rr_file_by_confidence{
	my $file_rr = shift;
	confess "No File RR!" if not -f $file_rr;
	my $seq = undef;
	if(rr_has_seq($file_rr)){
		$seq = seq_rr($file_rr);
	}
	my ($rr_sorted_fh, $rr_sorted_file) = tempfile();
	system_cmd("rm -f $rr_sorted_file");
	# Keep contact rows only
	system_cmd("sed -i '/^[A-Z]/d' $file_rr");
	system_cmd("sed -i '/^-]/d' $file_rr");
	# Some files have leading white spaces
	system_cmd("sed -i 's/^ *//' $file_rr");
	# Stable sort with -s option, i.e. maintain order in case confidence are equal
	# Also using -g option instead of -n because some predictors have exponential values in confidence
	system_cmd("sort -gr -s -k5 $file_rr > $rr_sorted_file");
	system_cmd("rm -f $file_rr");
	system_cmd("cat $rr_sorted_file >> $file_rr");
	system_cmd("rm -f $rr_sorted_file");
}

sub rr_has_seq{
	my $file_rr = shift;
	confess "ERROR! Input file $file_rr does not exist!" if not -f $file_rr;
	my $seq;
	open RR, $file_rr or confess "ERROR! Could not open $file_rr! $!";
	while(<RR>){
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		$_ =~ s/^\s+//;
		next if ($_ =~ /^PFRMAT/);
		next if ($_ =~ /^TARGET/);
		next if ($_ =~ /^AUTHOR/);
		next if ($_ =~ /^SCORE/);
		next if ($_ =~ /^REMARK/);
		next if ($_ =~ /^METHOD/);
		next if ($_ =~ /^MODEL/);
		next if ($_ =~ /^PARENT/);
		last if ($_ =~ /^TER/);
		last if ($_ =~ /^END/);
		last if ($_ =~ /^[0-9]/);
		$seq .= $_;
	}
	close RR;
	return 0 if not defined $seq;
	return 1;
}

sub xyz_pdb{
	my $chain = shift;
	confess "\nERROR! file $chain does not exist!" if not -f $chain;
	my %xyz_pdb = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$xyz_pdb{"".parse_pdb_row($_,"rnum")." ".parse_pdb_row($_,"aname")} = "".parse_pdb_row($_,"x")." ".parse_pdb_row($_,"y")." ".parse_pdb_row($_,"z");
	}
	close CHAIN;
	confess "\nERROR!: xyz_pdb is empty\n" if (not scalar keys %xyz_pdb);
	return %xyz_pdb;
}

sub system_cmd{
	my $command = shift;
	confess "EXECUTE [$command]" if (length($command) < 5  and $command =~ m/^rm/);
	system($command);
	if($? != 0){
		my $exit_code  = $? >> 8;
		confess "Failed executing [$command]!<br>ERROR: $!";
	}
}

sub seq_chain_with_gaps{
	my $chain = shift;
	my $flag = shift; # flag 1 if the left side dashes of the sequence are not wanted
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my $start = 1;
	# if flagged, keep start for trimming
	if (defined $flag){
		open CHAIN, $chain or confess $!;
		while(<CHAIN>){
			next if $_ !~ m/^ATOM/;
			if (parse_pdb_row($_,"rname") eq "GLY"){
				next if parse_pdb_row($_,"aname") ne "CA";
			}
			else{
				next if parse_pdb_row($_,"aname") ne "CB";
			}
			confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
			$start = parse_pdb_row($_,"rnum");
			last;
		}
		close CHAIN;
	}
	# 1.find end residue number
	my $end;
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		if (parse_pdb_row($_,"rname") eq "GLY"){
			next if parse_pdb_row($_,"aname") ne "CA";
		}
		else{
			next if parse_pdb_row($_,"aname") ne "CB";
		}
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		$end = parse_pdb_row($_,"rnum");
	}
	close CHAIN;
	# 2.initialize
	my $seq = "";
	for (my $i = 1; $i <= $end; $i++){
		$seq .= "-";
	}
	# 3.replace with residues
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		if (parse_pdb_row($_,"rname") eq "GLY"){
			next if parse_pdb_row($_,"aname") ne "CA";
		}
		else{
			next if parse_pdb_row($_,"aname") ne "CB";
		}
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		my $rnum = parse_pdb_row($_,"rnum");
		$rnum =~ s/[A-G]//g;
		substr $seq, ($rnum - 1), 1, $AA3TO1{parse_pdb_row($_,"rname")};
	}
	close CHAIN;
	confess "$chain has less than 1 residue!" if (length($seq) < 1);
	return (substr $seq, $start - 1);
}

sub seq_chain{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my $seq = "";
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		if (parse_pdb_row($_,"rname") eq "GLY"){
			next if parse_pdb_row($_,"aname") ne "CA";
		}
		else{
			next if parse_pdb_row($_,"aname") ne "CB";
		}
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		my $res = $AA3TO1{parse_pdb_row($_,"rname")};
		$seq .= $res;
	}
	close CHAIN;
	confess "$chain has less than 1 residue!" if (length($seq) < 1);
	return $seq;
}

sub parse_pdb_row{
	my $row = shift;
	my $param = shift;
	my $result;
	$result = substr($row,6,5) if ($param eq "anum");
	$result = substr($row,12,4) if ($param eq "aname");
	$result = substr($row,16,1) if ($param eq "altloc");
	$result = substr($row,17,3) if ($param eq "rname");
	$result = substr($row,22,5) if ($param eq "rnum");
	$result = substr($row,26,1) if ($param eq "insertion");
	$result = substr($row,21,1) if ($param eq "chain");
	$result = substr($row,30,8) if ($param eq "x");
	$result = substr($row,38,8) if ($param eq "y");
	$result = substr($row,46,8) if ($param eq "z");
	confess "Invalid row[$row] or parameter[$param]" if (not defined $result);
	$result =~ s/\s+//g;
	return $result;
}

sub res_num_res_name{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %rnum_rname = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$rnum_rname{parse_pdb_row($_,"rnum")} = parse_pdb_row($_,"rname");
	}
	close CHAIN;
	confess ":(" if not scalar keys %rnum_rname;
	return %rnum_rname;
}

sub calc_dist{
	my $x1y1z1 = shift;
	my $x2y2z2 = shift;
	my @row1 = split(/\s+/, $x1y1z1);
	my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
	my @row2 = split(/\s+/, $x2y2z2);
	my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
	my $d = sprintf "%.3f", sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
	return $d;
}

sub return_cb_or_ca_atom{
	my $seq = shift;
	my $rnum = shift;
	my $residue = substr $seq, $rnum - 1, 1;
	confess "rnum not defined!" if not defined $rnum;
	confess "Could not find residue name for $rnum!" if not $residue;
	return "CA" if $residue eq "G";
	return "CB";
}

sub seq_rr{
	my $file_rr = shift;
	confess "ERROR! Input file $file_rr does not exist!" if not -f $file_rr;
	my $seq;
	open RR, $file_rr or confess "ERROR! Could not open $file_rr! $!";
	while(<RR>){
		chomp $_;
		$_ =~ s/\r//g; # chomp does not remove \r
		$_ =~ s/^\s+//;
		$_ =~ s/\s+//g;
		next if ($_ =~ /^>/);
		next if ($_ =~ /^PFRMAT/);
		next if ($_ =~ /^TARGET/);
		next if ($_ =~ /^AUTHOR/);
		next if ($_ =~ /^SCORE/);
		next if ($_ =~ /^REMARK/);
		next if ($_ =~ /^METHOD/);
		next if ($_ =~ /^MODEL/);
		next if ($_ =~ /^PARENT/);
		last if ($_ =~ /^TER/);
		last if ($_ =~ /^END/);
		# Now, I can directly merge to RR files with sequences on top
		last if ($_ =~ /^[0-9]/);
		$seq .= $_;
	}
	close RR;
	confess "Input RR file does not have sequence row!</br>RR-file : <b>$file_rr</b></br>Please make sure that all input RR files have sequence headers!" if not defined $seq;
	return $seq;
}

sub print_usage{
	my $param_info = <<EOF;
--------------------------------------------
CONEVA v1.0 - Contact Assessment Toolkit
--------------------------------------------

PARAM              DEFAULT     DESCRIPTION
pdb   : file     : -         : Native PDB structure (mandatory)
rr    : file/dir : -         : RR file or directory
atom  : string   : cb        : Atom type (ca/cb/all/heavyatoms)
d     : float    : 8.0       : Distance threshold
smin  : int      : 24        : Minimum sequence separation (6 for short-range and 12 for medium-range)
smax  : int      : INF       : Maximum sequence separation (11 for short-range and 23 for medium-range)
jacc  : string   : N/A       : Contact selection for Jaccard Similarity (5, L/10, L/5, L/2, L, 2L)
all   : flag     : 0         : Show stats for all input contacts (slower)

Example:
./coneva.pl -rr ./1a3aA.rr -pdb 1a3aA.pdb
EOF
	print $param_info;
	print "\n";
	exit 1;
}
