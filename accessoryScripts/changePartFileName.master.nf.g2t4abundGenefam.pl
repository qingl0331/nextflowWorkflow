#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename qw(basename);
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage= "perl changePartFileName.g2t4abundGenefam.pl filelDir SampleSequencingDone_Master_190409missing15411X45.txt\nthis is to change the metagenome ID to metatranscriptome for matched patients' abundance and genefamilies.tsv\n";
my $inputDir = $ARGV[0];
my $keyFile = $ARGV[1];
die $usage unless $keyFile;
my%name;
open FH,$keyFile;
while(defined (my$line=<FH>)){
        next if $line=~ /^ShortID/;
        next if $line=~ /15412/;
        my@idv=split /\t+/,$line;
	next if $idv[4] eq "NA";
	next if ($idv[3] eq "NA" or $idv[3] eq "FAIL");
        $name{$idv[3]}=$idv[4];
}
close FH;

my@files=<$inputDir/*tsv>;
foreach my $file  (@files){
	my$fileName=basename $file;
	my($sName)=$fileName=~ /(1\d+X\d+)_\S+/;
	my($type)=$fileName=~ /L\d+_(\S+)\.tsv/;
	next if $type eq "pathcoverage";
	if(defined $name{$sName}){
		system("cp $file $name{$sName}.norm.$type.tsv");
	}
}
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------





