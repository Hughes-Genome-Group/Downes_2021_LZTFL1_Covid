#!/bin/bash -l
#$ -cwd

### (C) Damien Downes 2020.
### (C) Jelena Telenius 20th Nov 2018 (portability support)

# If you want more stringency, you can uncomment the below line
# To crash whenever anything within this run breaks ..

# set -e -o pipefail

# This is beneficial, as this level is just a wrapper.
# If any of the underlying perls DIE, this script will now abort at that point
# This may be a bit too stringent for certain points of the BASH parts,
# as it will crash the whole run, if any BASH commands return non-zero exit status (if you ls unexistent files for example)
# Thus : to temporarily turn this off (to for example, ls unexistent files) , you can (when you need it) :
# set +e +o
# and after you have passed the code part where you "know you sometimes get non-zero exit statuses but it is harmless",
# you can turn it back on by:
# set -e -o pipefail

##########################################
###     Run start - welcome message    ###
##########################################

echo
echo "This is $0 by Damien Downes (C) 2020"
echo
echo "Running in machine "
hostname
echo "Run started at "
date
echo
echo "Input parameters: $@"
echo

##########################################
###     Where to fetch things from     ###
##########################################



# Get the place "where we are now" - i.e. where the script cis3way.sh was, when it was called.
topdir=$( dirname $( which $0 ) )
# shell script arguments always contain $0 . It is the place the script was called from.
# so it will contain /wherever/it/was/cis3way.sh
# It can be relative path (such as ../cis3way.sh) - so we need to be cunning here, to get the absolute location.
# which $0 returns the absolute path
# the $() construct means "output of the command within the brackets" - the part which goes to STDOUT not STDERR
# so, $( which $0 ) will return the "output of the command which $0"
# dirname $( which $0 ) returns the directory part of this.
# these commands are fully portable (unlike 'fp' command - which is only available in CBRG cluster)

scriptdir="${topdir}/bin"
confFolder="${topdir}/conf"


##########################################
###     Loading the environment        ###
##########################################

# Fetch the genome sizes files, etc.

. ${confFolder}/genomeBuildSetup.sh
. ${confFolder}/loadNeededTools.sh
. ${confFolder}/serverAddressAndPublicDiskSetup.sh
# The above calls all the stuff in the file, so it is natively available in this script as well.
# This includes all subroutines and all variables and their default values (but the above only contains subroutines)

# Load all toolkits
setPathsForPipe 

# Set all genome reference files
. ${scriptdir}/genomeSetters.sh
CaptureDigestPath="UNDEFINED"
supportedGenomes=()
UCSC=()
setGenomeLocations

echo 
echo "Supported genomes : "
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do echo -n "${supportedGenomes[$g]} "; done
echo 
echo

echo "Reading public data area and server location for data hubs .."

SERVERTYPE="UNDEFINED"
SERVERADDRESS="UNDEFINED"
REMOVEfromPUBLICFILEPATH="NOTHING"
ADDtoPUBLICFILEPATH="NOTHING"
tobeREPLACEDinPUBLICFILEPATH="NOTHING"
REPLACEwithThisInPUBLICFILEPATH="NOTHING"

. ${confFolder}/serverAddressAndPublicDiskSetup.sh

setPublicLocations

echo
echo "SERVERTYPE ${SERVERTYPE}"
echo "SERVERADDRESS ${SERVERADDRESS}"
echo "ADDtoPUBLICFILEPATH ${ADDtoPUBLICFILEPATH}"
echo "REMOVEfromPUBLICFILEPATH ${REMOVEfromPUBLICFILEPATH}"
echo "tobeREPLACEDinPUBLICFILEPATH ${tobeREPLACEDinPUBLICFILEPATH}"
echo "REPLACEwithThisInPUBLICFILEPATH ${REPLACEwithThisInPUBLICFILEPATH}"
echo


#################################
###     Helper subs           ###
#################################

. "${scriptdir}/helpersubs.sh"
. "${scriptdir}/inputoutputtesters.sh"
# The above calls all the stuff in the file, so it is natively available in this script as well.
# This includes all subroutines and all variables and their default values (but the above only contains subroutines)

#################################
###     Default Parameters    ###
#################################

REenzyme="dpnII"
analysis="3way"

###############################
###     Input Parameters    ###
###############################


OPTS=`getopt -o g:: --long analysis:,path:,sample1:,sample2:,sample3:,sample4:,sample5:,sample6:,directories:,name:,genome:,restrictionenzyme:,version:,parameters:,annotation:,public:,frags: --  "$@"`

eval set -- "$OPTS"

while true ; do
    case "$1" in
        --analysis) analysis=$2 ; shift 2 ;;
        --path) path=$2 ; shift 2 ;;
        --sample1) sample1=$2 ; shift 2 ;;
        --sample2) sample2=$2 ; shift 2 ;;
        --sample3) sample3=$2 ; shift 2 ;;
        --sample4) sample4=$2 ; shift 2 ;;
        --sample5) sample5=$2 ; shift 2 ;;
        --sample6) sample6=$2 ; shift 2 ;;
        --directories) directories=$2 ; shift 2 ;;
        --name) name=$2 ; shift 2 ;;
        -g|--genome) genome=$2 ; shift 2 ;;
        --restrictionenzyme) REenzyme=$2 ; shift 2 ;;
        --version) version=$2 ; shift 2 ;;
        --parameters) parameters=$2 ; shift 2 ;;
        --annotation) annotation=$2 ; shift 2 ;;
        --frags) CaptureDigestPath=$2 ; shift 2 ;;
        --public) PublicPath=$2 ; shift 2 ;;
        --) shift; break;;
    esac
done

echo
echo "Parsed parameters are..."

# If '--analysis' given, but empty value, defaulting to 3way
if [ "${analysis}" = "" ]
    then
    analysis="3way"
    fi

# However, in all situations where sample3 is missing, defaulting to 2way
if [ "${sample3}" == "" ]
    then
    analysis="2way"
    fi


echo "Analysis Mode: ${analysis}"

echo "Path to input directories: ${path}"
echo "Sample1: ${sample1}"
echo "Sample2: ${sample2}"
echo "Sample3: ${sample3}"
echo "Sample4: ${sample4}"
echo "Sample5: ${sample5}"
echo "Sample6: ${sample6}"
echo "Directory list: ${directories}"
echo "Run name: ${name}"
echo "Genome: ${genome}"
echo "CC Version: ${version}"
echo "RE Enzyme: ${REenzyme}"
echo "Path to viewpoint parameters file: ${parameters}"
echo "Path to hub location: ${PublicPath}"
echo "Path to RE Fragments: ${CaptureDigestPath}"

###############################################
###     Check that parameters are fine      ###
###############################################

testInputParameters
# This sub comes from inputoutputtesters.sh
# If any of the parameters is empty, contains ONLY whitespace, or has whitespace within the value (contains spaces etc) - crashes with easy-to-interpret error message.

testInputDirectories
# This sub comes from inputoutputtesters.sh
# Check that all the dirs exist (no typos in dir name and path, thus)
# Checks that the dirs can be found, in the required format
# ${path}/${name}/F6_greenGraphs_combined_${name}_${version}

#################################
###     Perl script paths     ###
#################################

### Full paths to the four perl run scripts.

normaliser="${scriptdir}/cis3way_nomaliser.pl"

union="${scriptdir}/cis3way_unionbdg.pl"

windower="${scriptdir}/cis3way_windower.pl"

plotScripter="${scriptdir}/cis3way_ggScripter.pl"

DEseq2Scripter="${scriptdir}/cis3way_DEseq2Scripter.pl"

DEseq2Parser="${scriptdir}/cis3way_DEseq2Parser.pl"

peakyPrep="${scriptdir}/cis3way_PeakYprep.pl"

REcutter="${scriptdir}/${REenzyme}cutGenome4.pl"

#############################################
###     Config files - check and load     ###
#############################################

genomesSizes=""
setUCSCgenomeSizes
# This is in genomeSetters.sh - it uses the configuration loaded via genomeBuildSetup.sh
# it also checks that the file exists (aborts with helpful error message if not set properly)

# Checking if we have existing REdigest file,
# if not, generating here.

    printThis="Checking if this required fragment file is available: ${CaptureDigestPath}/${genome}_${REenzyme}_Fragments.txt"
    printToLogFile


if [ -s "${CaptureDigestPath}/${genome}_${REenzyme}_Fragments.txt" ]
    then
    frags="${CaptureDigestPath}/${genome}_${REenzyme}_Fragments.txt"
    printThis="Fragment file exists, digest not needed."
    printToLogFile
    else

    # Running the digestion ..
    printThis="Fragment file doesn't exist, digest required \n Running whole genome fasta ${REenzyme} digestion on ${genome}."
    printToLogFile

    echo "Run command: perl ${REcutter} ${GenomeFasta}"
    
    GenomeFasta=""
    setGenomeFasta
    # This comes from bin/genomeSetters.sh which reads the config from config/genomeBuildSetup.sh

    perl ${REcutter} "${GenomeFasta}"

    testedFile="genome_${REenzyme}_Fragments.txt"

    doTempFileTesting
 
    mv ${testedFile} "${genome}_${REenzyme}_Fragments.txt"
    madefile="${genome}_${REenzyme}_Fragments.txt"
    frags="$(pwd)/${madefile}"
    printThis="Generated fragment file: ${frags}"
    printToLogFile    
    fi
 

#############################################
###     Server and public area setup      ###
#############################################

echo "Parsing the public data area and server locations .."

PublicPath="${PublicPath}/${name}/"

echo "Updated PublicPath (disk path) to be : ${PublicPath}"

# Here, parsing the data area location, to reach the public are address..
diskFolder=${PublicPath}
serverFolder=""   
echo
parsePublicLocations
echo

tempJamesUrl="${SERVERADDRESS}/${serverFolder}"
JamesUrl=$( echo ${tempJamesUrl} | sed 's/\/\//\//g' )
ServerAndPath="${SERVERTYPE}://${JamesUrl}"

# Check the parses - the public area ownership test comes later (when the dir actually gets generated)

checkThis="${PublicPath}"
checkedName='${PublicPath}'
checkParse

checkThis="${ServerAndPath}"
checkedName='${ServerAndPath}'
checkParse

# ----------------------------------------------

#########################################
###     Normalisation run command     ###
#########################################

printThis="Running normalisation for the CaptureC gff files"
printToLogFile
# The above prints to both error and output logs. The sub resides in helpersubs.sh.
# This is beneficial, so the errors have some knowledge "whereabouts they happened"

echo
echo "Run command: perl ${normaliser} -path ${path} -dir ${directories} -viewpoints ${parameters} -name ${name} -version ${version}"
echo

perl ${normaliser} -path ${path} -dir ${directories} -viewpoints ${parameters} -name ${name} -version ${version}


####################################
###     Sorting of bedgraphs     ###
####################################

printThis="Sorting bedgraphs .."
printToLogFile

echo
echo "Moving into directory with bedgraphs"
echo

cd ${name}_cis_analysis
cp ${parameters} Used_Viewpoints.txt                # Make a copy of the viewpoints file for use in CaptureSee

echo "Now working in:"
echo "$(pwd)"
echo "Generating sorted bedgraph files.."
echo

echo "Now sorting normalised file:"
for file in *normalised.bdg

    do
    
    echo

    echo $file
    
    sortedname=$( echo $file | sed 's/.bdg//' ) 
    
    sort -k1,1 -k2,2n $file > ${sortedname}_sorted_TEMP.bedgraph

done

echo "Now sorting raw file:"
for file in *raw.bdg

    do
    
    echo

    echo $file
    
    sortedname=$( echo $file | sed 's/_raw.bdg//' ) 
    
    sort -k1,1 -k2,2n $file > ${sortedname}_sorted_TEMP_raw.bedgraph

done

###########################################################
###     Make union bedgraphs commands then run them     ###
###########################################################

printThis="Running bedtools union bedgraph for all oligos"
printToLogFile

echo "Run command: perl ${union} -dir ${directories} -viewpoints ${parameters} -samples ${sample1},${sample2},${sample3},${sample4},${sample5},${sample6}"
echo

perl ${union} -dir ${directories} -viewpoints ${parameters} -samples ${sample1},${sample2},${sample3},${sample4},${sample5},${sample6}

# module load bedtools
# all tools are loaded centralised : via calling the setup in /config/loadNeededTools.sh

commands="bedtools_commands_TEMP.txt"

bash $commands

#############################################################
###     Generate a mean track and a subtraction track     ###
#############################################################

printThis="Generate a mean track and a subtraction track"
printToLogFile

echo "Now doing maths for..."
for file in *normalised.unionbdg

    do

    oligoname=$( echo $file | sed 's/_normalised.unionbdg//' )
    echo "${oligoname}...Calculating means..."

    #final samples is CD14 with 2 reps.
    awk '{printf "%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", $1,$2,$3,(($4+$5+$6)/3),(($7+$8+$9)/3),(($10+$11+$12)/3),(($13+$14+$15)/3),(($16+$17+$18)/3),(($19+$20)/2)}' $file > "${oligoname}_TEMP_means.unionbdg"
    
    echo "...Separating..."
    
    cut -f 1,2,3,4  ${oligoname}_TEMP_means.unionbdg > ${oligoname}_${sample1}_mean.bdg
    cut -f 1,2,3,5  ${oligoname}_TEMP_means.unionbdg > ${oligoname}_${sample2}_mean.bdg
    cut -f 1,2,3,6  ${oligoname}_TEMP_means.unionbdg > ${oligoname}_${sample3}_mean.bdg
    cut -f 1,2,3,7  ${oligoname}_TEMP_means.unionbdg > ${oligoname}_${sample4}_mean.bdg
    cut -f 1,2,3,8  ${oligoname}_TEMP_means.unionbdg > ${oligoname}_${sample5}_mean.bdg
    cut -f 1,2,3,9  ${oligoname}_TEMP_means.unionbdg > ${oligoname}_${sample6}_mean.bdg    
    
    
done
#Sample A minus Samples B

echo "...Subtracting..."
for file in *TEMP_means.unionbdg
    do
    oligoname=$( echo $file | sed 's/_TEMP_means.unionbdg//' )
    awk '{printf "%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", $1,$2,$3,$4-$5,$4-$6,$4-$7,$4-$8,$4-$9,$5-$6,$5-$7,$5-$8,$5-$9,$6-$7,$6-$8,$6-$9,$7-$8,$7-$9,$8-$9}' $file > "${oligoname}_TEMP_subs1.unionbdg"
    
    cut -f 1,2,3,4  ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample1}_min_${sample2}.bdg
    cut -f 1,2,3,5  ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample1}_min_${sample3}.bdg
    cut -f 1,2,3,6  ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample1}_min_${sample4}.bdg
    cut -f 1,2,3,7  ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample1}_min_${sample5}.bdg
    cut -f 1,2,3,8  ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample1}_min_${sample6}.bdg
    cut -f 1,2,3,9  ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample2}_min_${sample3}.bdg
    cut -f 1,2,3,10 ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample2}_min_${sample4}.bdg
    cut -f 1,2,3,11 ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample2}_min_${sample5}.bdg
    cut -f 1,2,3,12 ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample2}_min_${sample6}.bdg
    cut -f 1,2,3,13 ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample3}_min_${sample4}.bdg
    cut -f 1,2,3,14 ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample3}_min_${sample5}.bdg
    cut -f 1,2,3,15 ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample3}_min_${sample6}.bdg
    cut -f 1,2,3,16 ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample4}_min_${sample5}.bdg
    cut -f 1,2,3,17 ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample4}_min_${sample6}.bdg
    cut -f 1,2,3,18 ${oligoname}_TEMP_subs1.unionbdg > ${oligoname}_${sample5}_min_${sample6}.bdg

    done    
    
for file in *TEMP_means.unionbdg
    do  
    oligoname=$( echo $file | sed 's/_TEMP_means.unionbdg//' )
    awk '{printf "%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", $1,$2,$3,$5-$4,$6-$4,$7-$4,$8-$4,$9-$4,$6-$5,$7-$5,$8-$5,$9-$5,$7-$6,$8-$6,$9-$6,$8-$7,$9-$7,$9-$8}' $file > "${oligoname}_TEMP_subs2.unionbdg"
    
    cut -f 1,2,3,4  ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample2}_min_${sample1}.bdg
    cut -f 1,2,3,5  ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample3}_min_${sample1}.bdg
    cut -f 1,2,3,6  ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample4}_min_${sample1}.bdg
    cut -f 1,2,3,7  ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample5}_min_${sample1}.bdg
    cut -f 1,2,3,8  ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample6}_min_${sample1}.bdg
    cut -f 1,2,3,9  ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample3}_min_${sample2}.bdg
    cut -f 1,2,3,10 ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample4}_min_${sample2}.bdg
    cut -f 1,2,3,11 ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample5}_min_${sample2}.bdg
    cut -f 1,2,3,12 ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample6}_min_${sample2}.bdg
    cut -f 1,2,3,13 ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample4}_min_${sample3}.bdg
    cut -f 1,2,3,14 ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample5}_min_${sample3}.bdg
    cut -f 1,2,3,15 ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample6}_min_${sample3}.bdg
    cut -f 1,2,3,16 ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample5}_min_${sample4}.bdg
    cut -f 1,2,3,17 ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample6}_min_${sample4}.bdg
    cut -f 1,2,3,18 ${oligoname}_TEMP_subs2.unionbdg > ${oligoname}_${sample6}_min_${sample5}.bdg

    done    

echo "...done."
echo       

#################################################
###     Convert all bredgraphs to bigwigs     ###
#################################################

printThis="Convert bedgraphs to bigwigs"
printToLogFile

# module load ucsctools
# all tools are loaded centralised : via calling the setup in /config/loadNeededTools.sh

echo "Starting final conversion of..."
for file in *.bdg
    do


    echo "$file"

    bigwigname=$( echo $file | sed 's/.bdg//' )

    sort -k1,1 -k2,2n $file > ${bigwigname}_TEMP_sorted.bdg

    echo
    cat ${bigwigname}_TEMP_sorted.bdg | grep -v chrUn > tmp
    bedGraphToBigWig tmp ${genomesSizes} ${bigwigname}.bw

done

echo "Removing all TEMP files"

rm -rf *TEMP*

#########################
###     Tidy up 1     ###
#########################

printThis="Cleaning up after ourselves .."
printToLogFile

mkdir 1_reports
mv *cisReport.txt* 1_reports/

mkdir 3_tracks

mkdir 3_tracks/A_replicates_raw
rm -rf *raw.bdg
mv *raw.bw 3_tracks/A_replicates_raw/

mkdir 3_tracks/B_replicates_normalised
rm -rf *normalised.bdg
mv *normalised*bw 3_tracks/B_replicates_normalised

mkdir 3_tracks/C_means
rm -rf *mean.bdg
mv *_mean*bw 3_tracks/C_means/

mkdir 3_tracks/D_subtractions
rm -rf *_min_*bdg
mv *_min_*bw 3_tracks/D_subtractions/


#########################################
###    Windowing for each viewpoint   ###
#########################################

printThis="Preparing binned and windowed tab files for use in R"
printToLogFile

    echo "Run command: perl ${windower} -viewpoints ${parameters} -samples ${sample1},${sample2},${sample3},${sample4},${sample5},${sample6}"
    echo

    perl ${windower} -viewpoints ${parameters} -samples ${sample1},${sample2},${sample3},${sample4},${sample5},${sample6}

#################################################
###    Generate Plotting R Scripts and run    ###
#################################################

printThis="Preparing R scripts for plotting"
printToLogFile

    echo "Run command: perl ${plotScripter} -viewpoints ${parameters} -annotation ${annotation}"
    echo
    
    perl ${plotScripter} -viewpoints ${parameters} -annotation ${annotation}

    echo "Running plotting script for..."

    for file in *plotting.R
        do
        viewpoint=$( echo $file | sed 's/_plotting.R//' )
        echo "${viewpoint}"
    
        Rscript $file
        
         
    done

#########################
###     Tidy up 2     ###
#########################

printThis="Cleaning up after ourselves .."
printToLogFile

mkdir 2_unionBedgraphs
mkdir 2_unionBedgraphs/B_normalised_counts
mv *normalised.unionbdg 2_unionBedgraphs/B_normalised_counts

mv  Sample_order.txt 2_unionBedgraphs/Sample_order.txt

mkdir 4_plotting

mkdir 4_plotting/A_parameters
mv *parameters.txt 4_plotting/A_parameters/

mkdir 4_plotting/B_binned
mv *bin.tab 4_plotting/B_binned/

mkdir 4_plotting/C_windowed
mv *window.tab 4_plotting/C_windowed/

mkdir 4_plotting/D_Rscripts
mv *.R 4_plotting/D_Rscripts/

mkdir 4_plotting/E_pdfs
mv *.pdf 4_plotting/E_pdfs/

##########################################
###     Make & Run DEseq2 Rscripts     ###
##########################################

printThis="Preparing R scripts for DEseq2"
printToLogFile

# Here changing to older version of R, to support the DESeq2 code part

changeRfrom3_4to3_2

    echo "Run command: perl ${DEseq2Scripter} --viewpoints ${parameters} -samples ${sample1},${sample2},${sample3},${sample4},${sample5},${sample6} -name ${name}"
    echo
    
    perl ${DEseq2Scripter} --viewpoints ${parameters} -samples ${sample1},${sample2},${sample3},${sample4},${sample5},${sample6} -name ${name}

printThis="Running DEseq2"
printToLogFile

    echo "Running DEseq2 script for..."

    for file in *.R
        do
        viewpoint=$( echo $file | sed 's/_DESeq2.R//' )
        echo "${viewpoint}"
    
        Rscript $file
           
    done

printThis="RParsing padj"
printToLogFile

    echo
    echo "RParsing padj for..."
    for file in *DESeq2_output.tab
        do
        analysis=$( echo $file | sed 's/_DESeq2_output.tab//' )
        echo "${analysis}"
    
        cut -f 1,7 $file > ${analysis}_TEMP_padj.tab
        
        perl ${DEseq2Parser}  --file ${analysis}_TEMP_padj.tab  --name ${analysis}
           
    done

printThis="Bedgraph to bigwig conversion"
printToLogFile

    echo
    echo "Starting final conversion of..."
    for file in *.bdg
        do
    
        echo "$file"
    
        bigwigname=$( echo $file | sed 's/.bdg//' )
    
        sort -k1,1 -k2,2n $file > ${bigwigname}_TEMP_sorted.bdg
    
        echo
        cat ${bigwigname}_TEMP_sorted.bdg | grep -v chrUn > tmp
        bedGraphToBigWig tmp ${genomesSizes} ${bigwigname}.bw
    
    done
    
    echo
    echo "Removing all TEMP files"
    
    rm -rf *TEMP*
    rm -rf tmp
    
#########################
###     Tidy up 3     ###
#########################

printThis="Cleaning up after ourselves .."
printToLogFile

mkdir 5_DESeq2
mkdir 5_DESeq2/A_Rscripts
mv *_DESeq2.R 5_DESeq2/A_Rscripts/
mkdir 5_DESeq2/B_columnfile
mv *DESeq2.tsv 5_DESeq2/B_columnfile/
mkdir 5_DESeq2/C_inputMatricies
mv *DESeq2_in.tsv 5_DESeq2/C_inputMatricies/
mkdir 5_DESeq2/D_raw_output
mv *DESeq2_output.tab 5_DESeq2/D_raw_output/
mkdir 3_tracks/E_pvalues
rm -rf *logpadj.bdg
mv *.bw 3_tracks/E_pvalues/

#######################################
###     Prepare files for PeakY     ###             Script does nto wort yet. Fails to assigne correct Prey ID.
#######################################

printThis="Preparing raw counts for peakY"
printToLogFile

echo "Run command: perl ${peakyPrep} --viewpoints ${parameters} --samples ${sample1},${sample2},${sample3},${sample4},${sample5},${sample6} --name ${name} --genome ${genome} --REfragments ${frags} --enzyme ${REenzyme}"
echo

perl ${peakyPrep} --viewpoints ${parameters} --samples ${sample1},${sample2},${sample3},${sample4},${sample5},${sample6} --name ${name} --genome ${genome} --REfragments ${frags} --enzyme ${REenzyme}

mkdir 2_unionBedgraphs/A_raw_counts
mv *raw.unionbdg 2_unionBedgraphs/A_raw_counts
mkdir 6_PeakyInputs
mv *tsv 6_PeakyInputs/
mv *key.bed 6_PeakyInputs/


###################################
###     A hub would be nice     ###
###################################

printThis="Making a data hub and CaptureSee Files"
printToLogFile

mkdir -p ${PublicPath}
rawrun="$(pwd)"

cd ${PublicPath}

ln -s ${rawrun}/Used_Viewpoints.txt
ln -s ${rawrun}/2_unionBedgraphs/Sample_order.txt
ln -s ${rawrun}/1_reports/*cisReport.txt .

mkdir ${genome}
cd ${genome}/
mkdir unionbedgraphs_normalised
cd unionbedgraphs_normalised
ln -s ${rawrun}/2_unionBedgraphs/B_normalised_counts/*.unionbdg .

cd ${PublicPath}

mkdir ${genome}/deseq_output
cd ${genome}/deseq_output
ln -s ${rawrun}/5_DESeq2/D_raw_output/*.tab .

cd ${PublicPath}


echo  > hub.txt

	echo hub ${name} >> hub.txt
    echo shortLabel ${name} >> hub.txt
    echo longLabel ${name} CaptureCompare >> hub.txt
	echo genomesFile genomes.txt >> hub.txt
    echo email damien.downes@ndcls.ox.ac.uk >> hub.txt
    
echo > genomes.txt 

	echo genome ${genome} >> genomes.txt
    echo trackDb ${genome}/tracks.txt >> genomes.txt
 
### Make a parent track of all the means

echo  > tracks.txt

	echo track Means >> tracks.txt
	echo type bigWig >> tracks.txt
	echo container multiWig >> tracks.txt
	echo shortLabel Mean >> tracks.txt
	echo longLabel CC_Means >> tracks.txt
	echo visibility hide >> tracks.txt
	echo aggregate transparentOverlay >> tracks.txt
	echo showSubtrackColorOnUi on >> tracks.txt
	echo maxHeightPixels 500:100:0 >> tracks.txt
	echo priority 1 >> tracks.txt
	echo   >> tracks.txt
	echo   >> tracks.txt

mkdir -p ${genome}/means

ln -s ${rawrun}/3_tracks/C_means/*.bw .

for file in *.bw
do
    
    trackname=$( echo $file | sed 's/.bw//' )

    echo track "${trackname}" >> tracks.txt
    echo bigDataUrl http://userweb.molbiol.ox.ac.uk"${PublicPath}"/"${genome}"/means/"${file}" >> tracks.txt
    echo shortLabel "${trackname}" >> tracks.txt
    echo longLabel "${trackname}" >> tracks.txt
    echo type bigWig >> tracks.txt
    echo parent Means >> tracks.txt
    echo color 0,0,0 >> tracks.txt
    echo   >> tracks.txt


    sed -i 's/\.\///' tracks.txt

done

mv *.bw ${genome}/means/


### Make a parent track of all the subtractions

	echo   >> tracks.txt
	echo   >> tracks.txt
	echo track Subtractions >> tracks.txt
	echo type bigWig >> tracks.txt
	echo container multiWig >> tracks.txt
	echo shortLabel Subs >> tracks.txt
	echo longLabel CC_Subs >> tracks.txt
	echo visibility hide >> tracks.txt
	echo aggregate transparentOverlay >> tracks.txt
	echo showSubtrackColorOnUi on >> tracks.txt
	echo maxHeightPixels 500:100:0 >> tracks.txt
	echo priority 2 >> tracks.txt
	echo   >> tracks.txt
	echo   >> tracks.txt

mkdir -p ${genome}/subtractions

ln -s ${rawrun}/3_tracks/D_subtractions/*.bw .

for file in *.bw
do
    
    trackname=$( echo $file | sed 's/.bw//' )

    echo track "${trackname}" >> tracks.txt
    echo bigDataUrl http://userweb.molbiol.ox.ac.uk"${PublicPath}"/"${genome}"/subtractions/"${file}" >> tracks.txt
    echo shortLabel "${trackname}" >> tracks.txt
    echo longLabel "${trackname}" >> tracks.txt
    echo type bigWig >> tracks.txt
    echo parent Subtractions >> tracks.txt
    echo color 0,0,0 >> tracks.txt
    echo   >> tracks.txt


    sed -i 's/\.\///' tracks.txt

done

mv *.bw ${genome}/subtractions/


### Make a parent track of all the pValues

	echo   >> tracks.txt
	echo   >> tracks.txt
	echo track pValues >> tracks.txt
	echo type bigWig >> tracks.txt
	echo container multiWig >> tracks.txt
    
	echo shortLabel pValues >> tracks.txt
	echo longLabel CC_pValues >> tracks.txt
	echo visibility hide >> tracks.txt
	echo aggregate transparentOverlay >> tracks.txt
	echo showSubtrackColorOnUi on >> tracks.txt
	echo maxHeightPixels 500:100:0 >> tracks.txt
	echo priority 3 >> tracks.txt
	echo   >> tracks.txt
	echo   >> tracks.txt

mkdir -p ${genome}/pvalues

ln -s ${rawrun}/3_tracks/E_pvalues/*.bw .

for file in *.bw
do
    
    trackname=$( echo $file | sed 's/.bw//' )

    echo track "${trackname}" >> tracks.txt
    echo bigDataUrl http://userweb.molbiol.ox.ac.uk"${PublicPath}"/"${genome}"/pvalues/"${file}" >> tracks.txt
    echo shortLabel "${trackname}" >> tracks.txt
    echo longLabel "${trackname}" >> tracks.txt
    echo type bigWig >> tracks.txt
    echo parent pValues >> tracks.txt
    echo color 0,0,0 >> tracks.txt
    echo   >> tracks.txt


    sed -i 's/\.\///' tracks.txt

done

mv *.bw ${genome}/pvalues/

mv tracks.txt ${genome}/

echo "Your hub is here:"
echo "${ServerAndPath}/hub.txt"


echo
echo "All done!"
echo

exit
