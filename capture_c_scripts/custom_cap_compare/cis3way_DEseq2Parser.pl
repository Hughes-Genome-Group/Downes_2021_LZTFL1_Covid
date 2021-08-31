#!/usr/bin/perl
use strict;
use Cwd;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;


### (C) Damien Downes 16th May 2018.

&GetOptions
(
    "file=s"=>\ my $input,   # --file
    "name=s"=>\ my $name,   # --file 
);

open (OUTPUT, ">$name\_logpadj.bdg") or die "cant open $name\_logpadj.bdg";
open (PADJ, $input) or die "Cannot open $input\n";
while (my $pvalue = <PADJ>)
	{
		chomp($pvalue);
		my ($frag, $padj) = split(' ', $pvalue);
        if ($frag =~ /\"baseMean\"/) {next;}
        if ($padj =~ /NA/) {next;}
        my $logpadg = -log($padj)/log(10);
        my $frag2 = substr $frag, 1, -1;
        my ($chr, $start, $stop) = split(/[:-]/,$frag2);
        if ($logpadg != 0)
        {
        print OUTPUT "$chr\t$start\t$stop\t$logpadg\n";
        }
    }
close PADJ;
close OUTPUT;
exit;
