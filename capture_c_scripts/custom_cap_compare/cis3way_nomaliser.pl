#!/usr/bin/perl -w
use strict;
use Cwd;
use Data::Dumper;
use Getopt::Long;
use Try::Tiny;


################################
##  Multi Run Cis Normaliser  ##
################################

# This script performs normalisation of all gff files listed in an oligo input for a set of samples and compares triplicates.
# Can function as a stand alone script or part of the CIS_analisis run shell which generates MEAN and SUBTRACTION tracks.

### (C) Damien Downes 16th May 2018 - a modified of Marieke Oudelaar's statistics script (capture_R_analysis1.pl) and cis analyser.


# Dependencies: 
#1. Input options for -path, -dir, -oligos, -name, -condition, -control -genome, -pf, -pu as described below
#2. Output folder will be a subdirectory of the folder that the script is in
#3. Input directories to be in this following structure:
#               .
#               |--Test_A
#               |  `--F6_greenGraphs_combined_Test_A_CC4
#               |     `COMBINED_CC4_OLIGO.gff
#               |
#               |--Test_B
#               |  `--F6_greenGraphs_combined_Test_B_CC4
#               |
#               |--Test_C
#               |  `--F6_greenGraphs_combined_Test_C_CC4
#               |
#               |--Control_A
#               |  `--F6_greenGraphs_combined_ControlA_CC4
#               |
#               |--Control_B
#               |  `--F6_greenGraphs_combined_ControlB_CC4
#               |
#               `--Control_C
#                  `--F6_greenGraphs_combined_ControlC_CC4


# Run Commands:

# Run with following options:       perl ./cis_normaliser.pl -path </full/path/to/CC4/output/> -dir test1,test2,test3,cont1,cont2,cont3 -oligo <oligo.file> -name <NAME>

# Example of a run command:         nohup perl cis_normaliser.pl -dir /t1-data1/WTSA_Dev/ddownes/Script_Testing/cis_normalisation/data -oligo oligos.txt -name run &
 
 
&GetOptions
(
    "path=s"=>\ my $path, 		        # -path         Path that contains all 6 Subdirectories (=the folders to all your different samples)
                                        #               generated by CC4 for the experiment you want to analyse
    "dir=s"=> \my @dirs,                # -dir          The directories that contain the gff files for all the samples (subdirectories) that you want to analyse,
                                        #               separated by only a comma (no space); these subdirectories need to contain the gff files, so it will be the ???F6???
                                        #               folder generated by CC4; you don???t have to give the full path to these folders, just direct down from the path
                                        #               you???ve entered (it will be two directories down)
    "viewpoints=s" => \my $vp_file,     # -viewpoints   VIEWPOINT	CHR VP_START VP_STOP EXCLSTART EXCLSTOP REGIONSTART REGIONSTOP BINSIZE WINDOWSIZE
                                        #               Note: chr is numeric  (1 NOT chr1)
    "name=s"=> \my $run_name,           # -name         Name of the folder containing the output files, can be whatever you like
    "version=s"=> \my $version,         # -version      Version of the analyser run.
);
 
 
@dirs = split(/,/,join(',',@dirs));
    
# Creates a folder for the output files - this will be a subdirectory of the file that the script is in
my $current_directory = cwd;
my $output_path= "$current_directory/$run_name\_cis_analysis";
if (-d $output_path){}
else {mkdir $output_path};

##### Doing analysis for each sample

my $name;

foreach my $name (@dirs)            #### Works through the array one at a time assigning the given directory to $name       
    {
    print "Current sample: $name \n";
    my $cis_summary_out = "$output_path\/$name\_cisReport.txt";
    open(VIEW_P, $vp_file) or die "Can't open $vp_file file";
    open (OUT_SUMMARY, ">$cis_summary_out") or die "can't open output file $cis_summary_out";
    print OUT_SUMMARY "Oligo\tCis\tTrans\tTotal\tPercent_Cis\n";
    
    while (my $target = <VIEW_P>)
        { 
        chomp $target;
        my ($viewID, $cis_chr, $vp_start, $vp_stop, @rest) = split(' ', $target);
        my $gff =  "$path/$name/F6_greenGraphs_combined_$name\_$version\/COMBINED_$version\_$viewID.gff";                  ## Will go balls up if folder path isn't exact!  - Common fail point                                                        
        open (FH, $gff) or die "can't open $name $viewID gff file which was supposed to be found in $gff , ";                                                             
        my $raw_out = "$output_path\/$name\_$viewID\_raw.bdg";               
        open(RAW, ">$raw_out") or die "Can't open $raw_out file";
        #### Calculating number of cis interactions, printing ALL raw counts to a bedgraph
                                                                                           
        my $cis_counter = 0;
        my $trans_counter = 0;
        while (my $line = <FH>)
                {
                chomp $line;
                my ($chr_test, $CC, $VP, $start2, $stop2, $value, $plus, $zero, $dot) = split(' ', $line);      
                print RAW "$chr_test\t$start2\t$stop2\t$value\n";               
                if ($chr_test eq $cis_chr)         #### FUTURE: Could add value in here for window to Mb region -- Could do with "if" comparing to co-ord ??frag start/stop. Add optional flag.
                        {
                        $cis_counter = $cis_counter + $value;
                        }
                else
                        {
                        $trans_counter = $trans_counter + $value;
                        }
                }
        my $fraction = sprintf("%.5f", ($cis_counter/($cis_counter+$trans_counter)));
        my $total = ($cis_counter+$trans_counter);
        print OUT_SUMMARY "$viewID\t$cis_counter\t$trans_counter\t$total\t$fraction\n";
        close FH;  
        $fraction = 0;
    
        #### Generating a normalised bdg of interactions per 100,000 unique cis interactions     Could potentially add super optional flag for different normalisation factor?
        
        my $cis_norm_out = "$output_path\/$name\_$viewID\_cis_normalised.bdg";               
        my $cis_denominator = $cis_counter / 100000;                                ## Will go balls up if cis = 0 - shouldn't happen unless input file error

        open(CISNORM, ">$cis_norm_out") or die "Can't open $cis_norm_out file";
        
        
        open (FH2, $gff) or die "can't re-open $viewID gff file";
        while (my $second_line = <FH2>)
                {
                chomp $second_line;
                my ($chr_test, $CC, $VP, $start2, $stop2, $value, $plus, $zero, $dot) = split(' ', $second_line);  
                if ($chr_test eq $cis_chr) ### Only prints reporters in cis - If adding Mb window, filter all counts outside region.   If we wanted to output trans could add flag option to generate that file.
                    {                      
                        my $norm_cis_value = $value / $cis_denominator;
                        print CISNORM "$chr_test\t$start2\t$stop2\t$norm_cis_value\n";
                    }
                }   
        close CISNORM;
        close FH2;                
        }
    
    close VIEW_P;
    close OUT_SUMMARY;

    }

exit;
