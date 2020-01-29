#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename qw(basename);
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage= "perl rf4lohImport.pl file \nthis is to reformat the loh output file for cBio portal import\n";
my $inputFile = $ARGV[0];
die $usage unless $inputFile;
my$fileName=basename $inputFile;
$fileName=~ s/\.targetWithGenes.txt//;
my (%hash,%header);
open my $fh,'+<',$inputFile;
while(defined (my$line=<$fh>)){
        my @array= split /\t/,$line;
        my$chr;
	if($array[0] eq "chrX"){
            $chr=23;
	}elsif($array[0] eq "chrY"){
            $chr=24;
        }else{
            ($chr)=$array[0]=~ /chr(\d+)/;
        }
        my$key="$chr.$array[1]";
        $hash{$key}="$array[0];$array[1]"."-"."$array[2];$array[6]";
        $header{$key}=$line;
}
close $fh;
open OUT1, ">$fileName.loh.import.txt";
open OUT2, ">$fileName.header.genes.txt";
print OUT2 "chr\tstart\tend\tname\tscore\tstrand\tgenes\n";
foreach my $ch(sort {$a<=>$b} keys %hash){
	print OUT1 $hash{$ch};
	print OUT2 $header{$ch};
}
close OUT1;
close OUT2;
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------





