#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename qw(basename);
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage= "perl changePartFileName.rfLog.pl file \nthis is to reformat sortmerna no_rRNA log file and change the file name\n";
my $inputFile = $ARGV[0];
die $usage unless $inputFile;
my$fileName=basename $inputFile;
$fileName=~ s/\.no_rRNA\.log//;
open OUT, ">$fileName.rf.log";
open my $fh,'+<',$inputFile;
while(defined (my$line=<$fh>)){
	if($line=~ /bowtie2 -q/){
		my ($id)=$line=~ /(1\w+)\.no_rRNA_humann2_temp/;
                $line=~ s/\/tmp\w+ -S/\/$id -S/;
	}
        print OUT $line;
}
close $fh;
close OUT;
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------





