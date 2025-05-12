#!/usr/bin/env perl
use v5.16.0;
use strict;
use warnings;

# EditDotplotic.pl rotate query_name +30, delete query_seq_sta, delete query_seq_end
# EditDotplotic.pl rotate query_name +30 Dotplotic.svg
# EditDotplotic.pl change ann_file.gff opacity 0.5

my $VERSION = "v1.0";

my @all_tags = qw(alignment scale_line1 scale_line2 query_name query_seq_sta query_seq_end sbjct_name sbjct_seq_sta sbjct_seq_end);

my $actions = {
	delete	=> {pn	=> [qw(target)],	target	=> [qw(query_name query_seq_sta query_seq_end sbjct_name sbjct_seq_sta sbjct_seq_end)],},
		# delete query_name
	move	=> {pn	=> [qw(target degree)],	target	=> [qw(query_name query_seq_sta query_seq_end sbjct_name sbjct_seq_sta sbjct_seq_end)],},
		# move query_name x+30
	rotate	=> {pn	=> [qw(target degree)],	target	=> [qw(query_name query_seq_sta query_seq_end sbjct_name sbjct_seq_sta sbjct_seq_end)],},
		# rotate query_name +30
	change	=> {pn	=> [qw(target attribute degree)],	target	=> [qw(alignment alignment_space scale_line1 scale_line2 query_name query_seq_sta query_seq_end sbjct_name sbjct_seq_sta sbjct_seq_end)],},
		# change alignment stroke-width 2
};

my @args = @ARGV;

my $config_file;
my $infile = "-";
if( ! $args[-1]){
	die "please input correctly!!\nSTOP\n";
} elsif(scalar(@args) == 1 && -f $args[0]){
    $config_file = $args[0];
} elsif(scalar(@args) == 2 && -f $args[0] && -f $args[1]){
    $config_file = $args[0];
	$infile = $args[1];
} elsif(-f $args[-1]){
	$infile = pop(@args);
}

my @commands;
if($config_file){
	@commands = read_in_config_file($config_file);
} else{
	@commands = read_commands(@args);
}

for my $command (@commands){
	version()		if($command eq "--version" || $command eq "-v");
	help($actions)	if($command eq "--help"    || $command eq "-h");
}

my $svglines = read_svgfile($infile);

my $outfmt;
for my $i (0..@commands-1){
	print STDERR $commands[$i]."\n";
	$svglines = one_command($svglines, $commands[$i]);
}

print $_ for @{$svglines};

exit;

#
#
#

sub read_svgfile{
	my $svgfile = shift;
	my $svg_lines;
	open my $IN, "<".$svgfile || die "Can't open $svgfile !!\nSTOP";
		push(@{$svg_lines}, <$IN>);
	close $IN;
	return $svg_lines;
}

sub read_commands{
	my @commands = split(/\,/, join(" ", @_));
	map{s/^\s+//g} @commands;
	return @commands;
}

sub one_command{
	my $innlines	= shift;
	my $command		= shift;

	my ($act, $target, $e1, $e2) = split(/\s/, $command);
    $target =~ s/\*/\(\\w\)/g;
	my ($degree, $attribute, $after); 
	if($act eq "delete"){
		;
	} elsif($act eq "move" || $act eq "rotate"){
		$degree = $e1;
	} elsif($act eq "change"){
		$attribute = $e1;
		$after = $e2;
	} else{
		print STDERR "\n\tYour action is ".$act.".\n";
		print STDERR "\tPlease choose an action from the following: ";
		print STDERR $_." " for sort keys %{$actions};
		print STDERR "\n\tSTOP\n";
		exit 0;
	}
	
	my %targted_tags;
	my $check;
	for my $tag (@{$actions->{$act}->{target}}){
		next if($target eq "alignment" && $tag eq "alignment_space");
		if($tag =~ /$target/){
			$targted_tags{$tag} = 1;
			$check = 1;
		}
	}
	die "\taction \'$act\' does not target \'$target\'!\n\tSTOP\n" if( ! $check);
	die "\taction \'$act\' does not target \'$target\'!\n\tSTOP\n" if($act ne "change" && ! scalar keys %targted_tags);

	my $outlines;
	for my $l (@{$innlines}){
		if($act eq "delete"){
			if($l =~ /\<\!\-\-(\w+)\-\-\>/){
				next if($targted_tags{$1});
			}
		} elsif($act eq "move" || $act eq "rotate"){
			if($l =~ /\<\!\-\-(\w+)\-\-\>/ && $targted_tags{$1}){
				my $tag = $1;
				if($act eq "move"){
					move($l, $degree);
				} elsif($act eq "rotate"){
					rotate($l, $degree, $tag);
				}
			}
		} elsif($act eq "change"){
			if($l =~ /\<\!\-\-(\S+)\-\-\>/ && ($targted_tags{$1} || $target eq $1)){
				my $tag = $1;
				change($l, $tag, $attribute, $after);
			}
		}

		push(@{$outlines}, $l);
	}

	return $outlines;
}

sub move{
	my $degree = $_[1];
	if($degree =~ /^([x|y])([\+\-])(\d+)$/){
		my ($xy, $sign, $val) = ($1,$2,$3); 
		if($_[0] =~ /$xy\=\"([\-\d\.]+)\"/){
			my $pre = $1;
			my $pst = ($sign eq "+")
				? $pre + $val
				: $pre - $val;
			$_[0] =~ s/$xy\=\"$pre\"/$xy\=\"$pst\"/;
		}
	}
}

sub rotate{ # $_[0] = $l;
	my ($degree, $tag) = ($_[1], $_[2]);
	if($degree =~ /^([\+\-])(\d+)$/){
		my ($sign, $val) = ($1,$2); 
		my ($r, $rx, $ry, $x, $y, $s);
		my $text = $1 if($_[0] =~ /\>(.*)\<\/text\>/);
		if($_[0] =~ /transform/){
			my $out = "\t<text ";
			while($_[0] =~ /([\-\w]+)\=\"([^\"]+)\"/){
				my ($a,$b) = ($1,$2);
				if($a eq "transform"){
					my $b2 = $b; $b2 =~ s/\s+//g;
					($r,$rx,$ry)=($1,$2,$3) if($b2 =~ /rotate\(([\-\d\.]+),([\-\d\.]+),([\-\d\.]+)/);
					die "Please check the SVG format!!\n\t".$_[0]."\nSTOP.\n" if(!$ry);
					my $r2 = ($sign eq "+")
						? $r + $val
						: $r - $val;
					$out .= "transform=\"rotate(".$r2.", ".$rx.", ".$ry.")\" ";
				} else{
					$out .= $a.'="'.$b.'" ';
				}
				$_[0] =~ s/$a\=//;
				$_[0]=~ s/$b//;
			}
			$_[0] = $out." >".$text."</text><!--".$tag."-->\n";
		} else{
			my $out = "\t<text ";
			my %contents;
			while($_[0] =~ /([\-\w]+)\=\"([^\"]+)\"/){
				my ($a,$b) = ($1,$2);
				$contents{$a} = $b;
				$out .= $a.'="'.$b.'" ';
				$_[0] =~ s/$a\=//;
				$_[0] =~ s/$b//;
			}
			$r = ($sign eq "+") ? $val : "-".$val ;
			$rx = $contents{"x"};
			$ry = $contents{"y"};
			$out .= "transform=\"rotate(".$r.", ".$rx.", ".$ry.")\" ";
			$_[0] = $out." >".$text."</text><!--".$tag."-->\n";
		}
	}
}

sub change{ # $_[0] = $l;
	my $tag		= $_[1];
	my $attribute	= $_[2];
	my $after	= $_[3];
	my $before = "";
	if($tag eq "alignment" && ($attribute ne "stroke-width" && $attribute ne "opacity") ){
		print STDERR "\n\tAttributes of the 'alignment' target are only 'stroke-width' or 'opacity'.\n\tSTOP.\n\n";
		exit;
	} elsif($tag eq "alignment_space" && ($attribute ne "fill") ){
		print STDERR "\n\tAttributes of the 'alignment_space' target are only 'fill'.\n\tSTOP.\n\n";
		exit;
	}

	if($attribute eq "fill" && is_valid_svg_color($after)){
		$after = "#".$after;
	}

	if($_[0] =~ /$attribute\=\"([\#\w\.\-]+)\"/){
		$before = $1;
		$_[0] =~ s/$attribute\=\"$before\"/$attribute\=\"$after\"/;
	} else{
		$_[0] =~ s/\>/$attribute\=\"$after\" \>/;
	}
}

sub is_valid_svg_color {
    my $color = shift;
    if ($color && $color =~ /^[0-9a-fA-F]{6}$/) {
        return 1;
    } else {
        return 0;
    }
}

sub read_in_config_file{
	my $file = shift;
	my @commands;
	open my $IN, "<".$file || die "Can't open the blast file!!\n";
		while(my $l = <$IN>){
			next if($l =~ /^#/ || $l !~ /\w/);
			chomp($l);
			$l = $1 if($l =~ /^([^\#]+)\#/);
			push(@commands, $l);
		}
	close $IN;
	return @commands;
}

sub help{
	my ($actions, $action_target);
	for my $act (sort keys %{$_[0]}){
		$actions .= $act.", " ; 
		my $target;
		$target = join(", ", @{$_[0]->{$act}->{target}});
		$action_target .= "\n\t".$act."\t:\t".$target;
	}
	chop($actions); chop($actions); 

	my $help_line="

	EditDotplotic.pl - a Perl program for modification of an SVG-format file from Dotplotic.pl


	The command of this program consists of an action, a target, an attribute, a change-value, and the input file:
	ex) \$ EditDotplotic.pl rotate query_name +30 Dotplotic.svg > Dotplotic.v2.svg

	This program can accept standard input and outputs to standard output:
	ex) \$ cat Dotplotic.svg | EditDotplotic.pl rotate query_name +30 > Dotplotic.v2.svg
	ex) \$ cat Dotplotic.svg | EditDotplotic.pl rotate query_name +30 | EditDotplotic.pl delete query_seq_sta > Dotplotic.v2.svg

	You can write multiple commands separated by a comma (,):
	ex) \$ EditDotplotic.pl rotate query_name +30 Dotplotic.svg, delete query_seq_sta, change sbjct_name font-size 20 > Dotplotic.v2.svg
	ex) \$ cat Dotplotic.svg | EditDotplotic.pl rotate query_name +30, delete query_seq_sta, change sbjct_name font-size 20 > Dotplotic.v2.svg


	Actions of this programs are as follows: $actions.

	Meanings of each targets are as follows:
	(query or sbjct)_name\t\t\t:\tname of query or subject
	(query or sbjct)_seq_(sta or end)\t:\tStart or End positions (base-pair) of query or subject
	scale_line(1 or 2)\t\t\t:\tScale lines

	Targets of each action are defined:$action_target

	The elemenets of the targets are as follows:
	<line> : alignment, scale_line1, scale_line2
	<text> : query_name, query_seq_sta, query_seq_end, sbjct_name, sbjct_seq_sta, sbjct_seq_end
	As for the details of attributes, such as 'stroke-width', 'fill', 'opacity' or 'font-size', please see about the SVG format.


	The 'change' action needs a target, an attribute of the target, and a value or word you choose:
	ex) \$ EditDotplotic.pl change query_name font-size 20 Dotplotic.svg > Dotplotic.v2.svg
	If you want to change attributes ('fill' and 'opacity') of annotaion, please write the name in the Dotplotic figure as the target:
	ex) \$ EditDotplotic.pl change Chr1_cds.gff opacity 0.5 Dotplotic.svg > Dotplotic.v2.svg
	CAUTION) attributes of the 'alignment' target have only 'stroke-width' or 'opacity';

	The 'delete' action only needs a target.
	ex) \$ EditDotplotic.pl delete query_name Dotplotic.svg > Dotplotic.v2.svg
	CAUTION) you can NOT get back targets once deleting:
	
	The 'move' action needs a target, and an offset. The degree consists of an axis('x' or 'y'), a direction ('+' or '-'), and a value:
	ex) \$ EditDotplotic.pl move query_name x-30 Dotplotic.svg > Dotplotic.v2.svg
	ex) \$ EditDotplotic.pl move sbjct_sta  y+15 Dotplotic.svg > Dotplotic.v2.svg

	The 'rotate' action needs a target, and an offseet. The degree consists of a direction ('+' or '-'), and a value:
	ex) \$ EditDotplotic.pl rotate query_name -30 Dotplotic.svg > Dotplotic.v2.svg

	";
	say $help_line;
	exit 1;
}

sub version{
	print STDERR "EditDotplotic $VERSION \n";
	exit 1;
}

