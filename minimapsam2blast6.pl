#!/usr/bin/env perl

=head1 NAME

minimapsam2blast6 - A program for convertion of SAM-format files from minimap2 into a Blast tabbular-format alignment file

=head1 SYNOPSIS

perl minimapsam2blast6 [-options] --sam <output SAM-format file of minimap2>

=head1 DESCRIPTION

This script converts an output file of blast search with "outfmt 6" into a SVG format file for visualization.

=over 4

=item B<--h(elp)>

Detailed help

=item B<--v(ersion)>

Detailed the version of the script

=item B<--verbose>

If specified, detailed information about the process will be displayed.

=item B<--sam>

SAM-format alignment file from minimap2: If no file is specified, this program accepts input from standard output. (e.g., $ cat minimap2.sam | minimap2sam2blast6.pl > minimap2.blast6.tsv) 

=item B<--identity>

Calculation method of Identity: You can choose one of 'BLAST' or 'sub'.
BLAST) Identities will be calculated with BLAST method 
(Identity(%) = 100 - (sum_of_'M' - (value_of_'NM' - sum_of_'I' - sum_of_'D'))/(sum_of_'M' + sum_of_'I' + sum_of_'D')x100), 
sub) Identities will be calculated based on substitution rate 
(Identity(%) = 100 - (value_of_'NM' - sum_of_'I' - sum_of_'D')/sum_of_'M'x100) (Default: BLAST)

=item B<--skip_multi>

If specified, only the first alignment is showed. (Default: Off)

=back

=head1 SEE ALSO


#

=head1 COPYRIGHT

#

=head1 AUTHOR

 Hideyuki, MIYAZAWA <hmiyazawa0209@gmail.com>

=cut

#
# Module Dependence
#

use v5.16.0;
use strict;
use warnings;
use feature 'say';
use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);
use Pod::Text;
use File::Copy;
use File::Spec;
use File::Path;
use Data::Dumper;
use Cwd;

my $CLASS = "minimap2sam2blast6";
my $VERSION = "1.0.0";

my ($help, $version, $verbose);
my $sam_file //= "-";
my ($print_cigar, $skip_multi_mapping);
my $idn = "Blast";

GetOptions(
    'help|h'            => \$help,                  # Support --help or -h
    'version|v'         => \$version,               # Support --version or -v
    'sam=s'             => \$sam_file,              # 
    'identity|idn=s'    => \$idn,                   # 
    'verbose'           => \$verbose,               # 
    'cigar'             => \$print_cigar,           # 
    'skip_multi'        => \$skip_multi_mapping,    # 
) || die "Error in command-line arguments. Use --help for usage.\n";

# Display help
my $this_file = $0;
help($help, $this_file);

# Display version
version($version, $VERSION);

if($verbose){
    if($idn eq "Blast"){
        print STDERR "# Identity is calculated using Blast method.\n";
        print STDERR "# \tIdentity = 100 - (sum_of_'M'-(value_of_'NM'-sum_of_'I'-sum_of_'D'))/(sum_of_'M'+sum_of_'I'+sum_of_'D')*100\n";
    } elsif($idn eq "sub"){
        print STDERR "# Identity is calculated with substitution rate.\n";
        print STDERR "# \tIdentity = 100 - (value_of_'NM'-sum_of_'I'-sum_of_'D')/sum_of_'M'*100\n";
    } else{
        print STDERR "Unknown Identity style: $idn\nSTOP\n"; exit;
    }
    print STDERR "# All CIGARs are outputted at the end of line.\n" if($print_cigar);

    print STDERR "# output blast format : \"6 std qlen slen\"\n";
    print STDERR "# skipping multi mappped alignment\n" if($skip_multi_mapping);
    print STDERR "# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, query length, subject length\n";
    print STDERR "# CAUTION) All Evalues are outputted as 0.0\n";
}

my $check_flag = check_flag();
my $direction  = direction();

die "I cant open the file, ".$sam_file." !!\n" if($sam_file ne "-" && !-e $sam_file);
my (%align, %slen, %qlen, %ctg, %ctg2ref, %refbed, @queries, $minimap2); # $align{$query} = qw(samlines);
open my $IN, "<".$sam_file || die "Can't open the SAM-format file!!\n"; my $i=0;
while(my $l = <$IN>){
	chomp($l); ++$i;
	my @d = split(/\s+/, $l);
	if($d[0] =~ /^\@SQ/){
		my $sname = $1 if($d[1] =~ /SN\:(\S+)$/);
		my $slen  = $1 if($d[2] =~ /LN\:(\d+)$/);
		if(!$sname || !$slen){
            die "Can't find name and length of subjects. Please check the file!\nSTOP\n"
        }
		$slen{$sname} = $slen ;
	} elsif($d[0] =~ /^\@HD/ || $d[0] =~ /^\@RG/){
		;
	} elsif($d[0] =~ /^\@PG/){
		;
	} else{
		my ($query,$flag,$sbjct,$ssta,$mapq,$cigar)=($d[0],$d[1],$d[2],$d[3],$d[4],$d[5]);
		if (! $cigar){
			print STDERR $d[0]."!!!\n";
			print STDERR "line_n= ".$i."\t".substr($l, 0, 100)."\t".scalar(@d)."\n";
			print STDERR $d[$_]."\n" for(0..@d-1);
			die "Can't find CIGAR from this line. Please check the file!\nSTOP\n"
		}
		next if(! $check_flag->($d[1])); # skip if the flag means 'read unmapped'
		my $drc = $direction->($d[1]); # 'fwd' or 'rev'
		$ctg2ref{$query} = $sbjct ;
		push(@queries, $query) if(!$queries[0] || $queries[-1] ne $query) ;
		my ($qsta,$qend,$qlen,$total_M,$del_num,$del_siz,$ins_siz,$ins_num) = cigar_param($cigar, $drc);
		# In minimap2 and Blast, Query and Subject are inverted (see cigars outputted on 'blast -outfmt 17' and 'minimap2 -a').
		# So, exchange them.
		# ($del_num, $ins_num) = ($ins_num, $del_num);
		# ($del_siz, $ins_siz) = ($ins_siz, $del_siz);

		if( ! $qlen{$query}){
			$qlen{$query} = $qlen;
		} else{
			next if($skip_multi_mapping);
			next if($qlen{$query} != $qlen);
		}
		$del_num //= "0";

		my $mlen = $total_M + $del_siz + $ins_siz ;
		my $send = $ssta + $total_M + $del_siz -1 ;

		my $dp_scr	= ($l =~ /AS\:i\:([\-\d]+)\s/)	? $1 : "0" ;	# AS │  i   │ DP alignment score
		my $err_siz	= ($l =~ /NM\:i\:(\d+)\s/)		? $1 : "0" ;	# NM │  i   │ Total number of mismatches and gaps in the alignment
		next if($dp_scr < 0);

		my $ident = "0" ;
		if($idn eq "Blast"){
			$ident = ($mlen - $ins_siz - $err_siz)/$mlen * 100;
		} elsif($idn eq "sub"){
			# substitution size = $err_siz - $ins_siz - $del_siz ;
			$ident = ($total_M - ($err_siz - $ins_siz - $del_siz) )/$total_M * 100;
		}
		my $eval = "0.0";
		my $gap_n = $ins_num + $del_num;

		# blast outfmt "6 std qlen slen" + cigar
		($ssta,$send)=($send,$ssta) if($drc eq "rev");
		if ( ! $slen{$sbjct}){
			die "The length of $sbjct is unknown!!\nPlease check the SAM-format file!!\nSTOP\n";
		}
		print $_."\t" for($query, $sbjct, $ident, $mlen, $err_siz, $gap_n, $qsta, $qend, $ssta, $send, $eval, $dp_scr, $qlen, $slen{$sbjct});
		print $cigar."\t" if($print_cigar);
		say "";
	}
}
close $IN;

exit;

#
#
#

exit;

sub help{
    if ($_[0]) {
        my $parser = Pod::Text->new();  # Create a Pod::Text parser
        $parser->parse_from_file($0);   # Parse POD from the current file
        exit;
    }
    return 1;
}

sub version{
    if ($_[0]) {
        print "Dotplotic version $_[1]\n";
        exit;
    }
}

sub direction{ # closure : check_flag( flag )
	my @bi  = qw(1 2 4 8 16 32 64 128 256 512 1024 2048);
	my @badb= qw(1 1 1 1  0  1  1   1   1   1    1    1);
	return sub{
		my $flag = shift(@_);
		return "fwd" if($flag==0);
		my $out = 1 ;
		for my $i (0..11){
			my $j = 11-$i;
			if(!$bi[$j]){
                die "The flag is weird!!\nPlease check the SAM-format file!!\nSTOP\n";
            }
			if($flag>=$bi[$j]){
				$flag -= $bi[$j];
				$out *= $badb[$j];
			}
		}
		return ($out) ? "fwd" : "rev" ;
	}
}

sub check_flag{ # closure : check_flag( flag )
	my @bi  = qw(1 2 4 8 16 32 64 128 256 512 1024 2048);
	my @badb= qw(1 1 0 1  1  0  1   1   1   1    1    1);
	return sub{
		my $flag = shift(@_);
		return 1 if($flag==0);
		my $out = 1 ;
		for my $i (0..11){
			my $j = 11-$i;
			if(!$bi[$j]){
                die "The flag is weird!!\nPlease check the SAM-format file!!\nSTOP\n";
            }
			if($flag>=$bi[$j]){
				$flag -= $bi[$j];
				$out *= $badb[$j];
			}
		}
		return $out;
	}
}

sub cigar_param{ # $cigar
	my $cigar = $_[0];
	my $drc = $_[1];
	my ($qsta, $qend, $qlen) = (0,0,0);
	if($drc eq "fwd"){
		$qsta = ($cigar =~ /^(\d+)[S|H]/) ?  $1+1 : 1 ;
	} else{
		$qsta = ($cigar =~ /(\d+)[S|H]$/) ?  $1+1 : 1 ;
	}
	$qlen = $qsta;
	if($drc eq "fwd"){
		$qlen += ($cigar =~ /(\d+)[S|H]$/) ?  $1 : 0 ;
	} else{
		$qlen += ($cigar =~ /^(\d+)[S|H]/) ?  $1 : 0 ;
	}
	my ($total_M,$del_num,$del_siz,$ins_siz,$ins_num) = (0,0,0,0,0) ;
	while($cigar =~ /^(\d+)(\D)/){
		my ($v,$t) = ($1,$2);
		my $out = $v.$t;
		$total_M += $v if($t eq "M");
		# In minimap2, Query and Subject are inverted (see cigars on 'blast -outfmt 17' and 'minimap2 -a').
		# So, exchange them.
		if($t eq "D"){
			++$del_num;
			$del_siz += $v;
		}
		if($t eq "I"){
			++$ins_num;
			$ins_siz += $v;
		}
		$cigar =~ s/^$out//;
	}
	$qlen += $total_M + $ins_siz - 1 ;
	$qend = $qsta + $total_M + $ins_siz -1 ;
	return ($qsta,$qend,$qlen,$total_M,$del_num,$del_siz,$ins_siz,$ins_num);
}
