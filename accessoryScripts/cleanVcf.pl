#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename qw(basename);
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage= "perl cleanVcf.pl file \nthis is to clean normal vcf file\n";
my $inputFile = $ARGV[0];
die $usage unless $inputFile;
my$fileName=basename $inputFile;
$fileName=~ s/\.g.vcf//;
open my $fh, '<', $inputFile;
open my $out,'>',"$fileName.NormalDNA.clean.vcf";
while(defined (my$line=<$fh>)){
                next if $line=~ /\t\*\t/;
                print $out $line;
}
close $out;
close $fh;

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------





