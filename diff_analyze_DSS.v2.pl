#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use Cwd;
use File::Basename;
my $help = "";
my $data_path = "";
my $format ="";
my $deletemechrX = "no";
my $samplelist = "";
my $comparegroup = "";
my $dss = "no";
my $annotate = "no";
my $outdir = "out";

my $option = GetOptions(
	"help|h|?"          =>\$help,
	"data_path|d=s"     =>\$data_path,
	"format|f=s"        =>\$format,
	"deletechrX|a:s"    =>\$deletemechrX,
	"samplelist|s=s"    =>\$samplelist,
	"comparegroup|c=s"  =>\$comparegroup,
	"diff|i:s"          =>\$dss,
	"annotate|a:s"      =>\$annotate,
	"outdir|d:s"        =>\$outdir,
);

if($help or $data_path eq ""){
	print "This is diff_analyze_DSS.pl used to detect differentially methylated loci and region from single nucleotide resolution sequencing data.\n";
	print "Usage:\n";
	print "\t perl $0 -data_path cov -format cov --deletechrX yes -samplelist sample.group.xls -comparegroup malignant--benign --diff yes --annotate no --outdir out\n";
	print "\t\t -data_path|d          inputdir of count or cov file path\n";
	print "\t\t -format|f             inputdata format(count|cov)\n";
	print "\t\t --deletechrX|a         whether delete CpG sites on chrX when convert cov to count(yes|no,default:no)\n";
	print "\t\t -samplelist|s         samplelist with sample and group information\n";
	print "\t\t -comparegroup|c       compare group,eg:tumor--normal,should be consistent with samplelist\n";
	print "\t\t --diff|i              whether do dss analyze or not (yes|no,default:no)\n";
	print "\t\t --annotate|a          whether do annotate analyze or not(yes|no,default:no)\n";
	print "\t\t --outdir|d            output dir(default:./out)\n";
	exit;
}

my $cwd;
if ($0 =~ m{^/}) {
	$cwd = dirname($0);
}else{
	$cwd = dirname(getcwd()."/$0");
} 
my $Rscript="/anchordx-opt/local/bin/Rscript";
$data_path = abs_path($data_path);
$outdir = abs_path($outdir);
`mkdir -p $outdir` unless -d $outdir;

#------Convert cov 2 count---------------------------------------#
if ($format eq "cov"){
	$data_path = cov2count($data_path,$outdir);	
}

#------Get the sample and compare group information--------------#
my ($treatment,$control) = (split('--',$comparegroup));
my ($treatment_samples,$control_samples) = readsamplelist($samplelist,$treatment,$control);

#------DSS diff analyse------------------------------------------#
if ($dss eq "yes"){
	`$Rscript $cwd/DSS.R $data_path $treatment_samples--$control_samples $treatment--$control $outdir`;
}
#------annotate--------------------------------------------------#
if($annotate eq "yes"){
	my @file = `ls $outdir/$treatment--$control/*DMLtest.xls  $outdir/$treatment--$control/*.dml.*[0123456789].xls`;
	anno("dml",@file);
	my @file2 = `ls $outdir/$treatment--$control/*.dmr.*[0123456789].xls`;
	anno("dmr",@file2);
}

#------sub function----------------------------------------------#
sub cov2count{
	my ($path,$outdir) = @_;
	my $count_path = "$outdir/count";
	`mkdir -p $count_path` unless -d $count_path;
	my @infile = `ls $path/*gz`;
	my %call;
	my $covstatfile = "$count_path/coverage.stat.xls";
	open CF,">$covstatfile" or die $!;
	print "coverage_cpg_count\tfile\n";
	foreach(@infile){
		my $infile = $_;
		chomp $infile;
		my $outname = `ls $infile|awk -F "/" '{print \$NF}' |awk -F "_R1" '{print \$1}'`;
		chomp $outname;
		my $out = "$count_path/$outname".".count.xls";
		my $stat=0;
		open IN,"gzip -cd $infile|" or die $!;
			while(<IN>){
				chomp;
				my ($chr,$pos,$meth,$nometh) = (split("\t",$_))[0,1,4,5];
				my $flag = "$chr"."-"."$pos";
				my $total = $meth + $nometh;
				$call{$flag} = "$chr\t$pos\t$total\t$meth";	
			}
		close IN;
#		open BD,"/home/longhui_deng/script/CpGcount/methylation_PanCan9921Sites_padded_25-rm_LC_CR_outlier_161215_padded_100.cpg.bed" or die $!;
		open BD,"/home/longhui_deng/script/CpGcount/batronman/methylation_PanCan9921Sites_padded_25-rm_LC_CR_outlier_161215_padded_100.cpg.annovar.bed" or die $!;
		open OUT,">$out" or die $!;
		print OUT "chr\tpos\tN\tX\n";
		while(<BD>){
			chomp;
			my ($chr,$start) = (split("\t",$_))[0,2];
			#****20170803edit:delete chrX*************#
			if ($deletemechrX eq "yes"){
				if ($chr ne "chrX"){
					my $flag = "$chr"."-"."$start";
					if (exists $call{$flag}){
						print OUT "$call{$flag}\n";
						$stat++;
					}else{
						print OUT "$chr\t$start\t0\t0\n";
					}
				}
			}else{
				my $flag = "$chr"."-"."$start";
				if (exists $call{$flag}){
					print OUT "$call{$flag}\n";
				}else{
					print OUT "$chr\t$start\t0\t0\n";
				}
			}
		}	
		close BD;
		close OUT;
		print CF "$stat\t$infile\n";
	}
	close CF;
	return $count_path;
}

sub readsamplelist{
	my ($samplelist,$treatment,$control) = @_; 
	my ($treatment_samples,$control_samples);
	open SS,"$samplelist" or die $!;
	while(<SS>){
		chomp $_;
		my ($sample,$group) = split("\t",$_);
		if ($group eq $treatment){
			$treatment_samples .= ($treatment_samples)? (",".$sample):($sample);
		}elsif($group eq $control){
			$control_samples .= ($control_samples)?(",".$sample):($sample);
		}
	}
	return($treatment_samples,$control_samples);
}

sub anno{
	my ($flag,@file) = @_;
	foreach (@file){
		chomp;
		my $tmp = `ls $_|sed 's/.xls//g'`;
		chomp $tmp;
		if($flag eq "dml"){
			`less $_ |sed 's/:/\t/g'|awk 'NR>1{print \$1"\t"\$2-1"\t"\$2"\t0\t0"}' >$tmp.bed`;
		}else{
			`less $_ |sed 's/:/\t/g'|awk 'NR>1{print \$1"\t"\$2"\t"\$3"\t0\t0"}' >$tmp.bed`;
		}
		`perl $cwd/annotate_variation.pl -buildver hg19 $tmp.bed /Scratch/Database/annovar_db/humandb`;
		`echo "anno.type	anno.gene	chr	pos" > $tmp.anno.gene.xls`;
		`less $tmp.bed.variant_function |cut -f1-4 >>$tmp.anno.gene.xls`;
		if($flag eq "dml"){
			`less $tmp.xls |cut -f1 > $tmp.1`;
			`less $tmp.xls |cut -f2- > $tmp.2`;
			`less $tmp.anno.gene.xls|cut -f1,2 >$tmp.gene`;
			`paste $tmp.1 $tmp.gene $tmp.2 > $tmp.anno.xls`;
		}else{
			`less $tmp.xls |cut -f1,2,3 > $tmp.1`;
			`less $tmp.xls |cut -f4- > $tmp.2`;	
			`less $tmp.anno.gene.xls|cut -f1,2 > $tmp.gene`;
			`paste $tmp.1 $tmp.gene $tmp.2 > $tmp.anno.xls`;
		}
		my %type;
		my $sum;
		open TT,"$tmp.anno.xls" or die $!;
		while(<TT>){
			chomp;
			if($.==1){next}
			my $type = (split("\t",$_))[($flag eq "dml")? 1:3];
			if(!exists $type{$type}){
				$type{$type} = 1;
			}else{
				$type{$type}++;
			}
			$sum++;
		}
		close TT;

		open ST,">$tmp.anno.stat.xls" or die $!;
		print ST "anno.gene.type\tcount\tpercent\n";
		print ST "all\t$sum\t100%\n";
		for my $type (sort{$type{$b} <=> $type{$a}} keys %type){
			my $pct = sprintf("%.2f",$type{$type}/$sum*100);	
			print ST "$type\t$type{$type}\t$pct%\n";
		}
		close ST;
		if($flag eq "dml"){
			`Rscript $cwd/pie.R $tmp.anno.stat.xls $tmp.anno.stat.pdf DML`;
		}else{
			`Rscript $cwd/pie.R $tmp.anno.stat.xls $tmp.anno.stat.pdf DMR`;
		}
		`convert -density 300 $tmp.anno.stat.pdf $tmp.anno.stat.png && rm $tmp.anno.stat.pdf`;
		`rm -fr $tmp.bed $tmp.bed.exonic_variant_function $tmp.bed.variant_function $tmp.bed.log $tmp.anno.gene.xls $tmp.1 $tmp.2 $tmp.gene`;
	}
}





__END__













