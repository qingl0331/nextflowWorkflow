#!/usr/bin/perl -w
use strict;
use warnings;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage= "perl calcuAf.pl block.norm.id.NormalDNA.clean.vcf id.TumorDNA_Hg38_final.vcf outName\n";
my $nVcf = $ARGV[0];
my $tVcf = $ARGV[1];
my $outName = $ARGV[2];
die $usage unless $outName;
my%normHets;#{chr.pos}=[cov,baf]
open my $fh,'<',$nVcf;
while(defined (my$line=<$fh>)){
	next if $line=~ /^#/;
	next unless $line=~ /0\/1/;
	chomp $line;
	my @bafs=split /\t/,$line;
	my($dp)=$bafs[7]=~ /DP=(\d+)\;/;
	my($baf)=$bafs[7]=~ /AF=(\d*\.?\d*)\;/;
	my$chrPos="$bafs[0].$bafs[1]";
	$normHets{$chrPos}=[$dp,$baf];
	#print $out "$bafs[0]\t$bafs[1]\t$dp\t$baf\n";
}
open my $out2, '>',"$outName.tumor.baf.txt";
print $out2 "chr\tposition\tcoverage\tbaf\n";
open my $fh2,'<',$tVcf;
while(defined (my$line2=<$fh2>)){
        next if $line2=~ /^#/;
        chomp $line2;
        my @bafs2=split /\t/,$line2;
	my$chrPos2="$bafs2[0].$bafs2[1]";
	next unless defined $normHets{$chrPos2};
        my($dp2)=$bafs2[7]=~ /DP=(\d+)\;/;
        my($baf2)=$bafs2[7]=~ /AF=(\d*\.?\d*)\;/;
        print $out2 "$bafs2[0]\t$bafs2[1]\t$dp2\t$baf2\n";
}
close $out2;
open my $fh3, '<',"$outName.tumor.baf.txt";
my%tum;#{chr.pos}=[cov,baf, count]
while(defined (my$line3=<$fh3>)){
        next if $line3=~ /position/;
	chomp $line3;
	my @bafs3=split /\t/,$line3;
	#$bafs3[0]=~ s/chr//;
	my$chrPos3="$bafs3[0].$bafs3[1]";
	my$count=1;
	if(!defined $tum{$chrPos3}){
		$tum{$chrPos3}=[$bafs3[2],$bafs3[3],$count];
	}else{
		$tum{$chrPos3}->[0]+=$bafs3[2];
		$tum{$chrPos3}->[1]+=$bafs3[3];
		$tum{$chrPos3}->[2]+=$count;
	}
}
open my $out3, '>',"$outName.tumor.ave.baf.txt";
print $out3 "chr\tposition\tcoverage\tbaf\n";
foreach my $cp (sort keys %tum){
	my$covAve=int($tum{$cp}->[0]/$tum{$cp}->[2]);
	my$bafAve=$tum{$cp}->[1]/$tum{$cp}->[2];
	my@chpos=split /\./, $cp;
	print $out3 "$chpos[0]\t$chpos[1]\t$covAve\t$bafAve\n";
}
close $out3;
open my $out, '>',"$outName.normal.het.baf.txt";
print $out "chr\tposition\tcoverage\tbaf\n";
foreach my $chpo (sort  keys %normHets){
	next unless defined $tum{$chpo};
	my@chpos=split /\./, $chpo;
	print $out "$chpos[0]\t$chpos[1]\t$normHets{$chpo}->[0]\t$normHets{$chpo}->[1]\n";
}
close $out;
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub cleanVcf{
	my$gVcf=shift;
	my$outName=shift;
	open my $fh, '<',$gVcf;
	open my $out,'>',"$outName.NormalDNA.clean.vcf";
	while(defined (my$line=<$fh>)){
		next if $line=~ /\t\*\t/;
		print $out $line;
	}
	close $out;
}



