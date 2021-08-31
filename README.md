# Downes (2021) Identification of LZTFL1 as a candidate effector gene at a COVID-19 risk locus.


**Linkage Analysis**

Linkage analysis for variants in linkage with lead variants used [LD Proxy](https://ldlink.nci.nih.gov/?tab=ldproxy) and [LD Matrix](https://ldlink.nci.nih.gov/?tab=ldmatrix)


**Micro RNA binding**

3' UTR effects of rs35624553 were determined using [TargetScan](http://www.targetscan.org/vert_71/), [miRdSNP](http://mirdsnp.ccr.buffalo.edu/browse-genes.php), and [MicroSNiPer](http://vm24141.virt.gwdg.de/services/microsniper/index.php)


**Splicing**

SpliceAI masked splicing predictions were downloaded from the github [repository](https://github.com/Illumina/SpliceAI).

Splicing quantitative trait loci (sQTLs) were identified in the [GTEx Database](https://gtexportal.org/home/) 

**Chromatin accessibility and transcription factor binding prediction**

[deepHaem](https://github.com/rschwess/deepHaem) was implemented using *run_deephaem_4k_predictions.sh*

[Sasquatch](https://apps.molbiol.ox.ac.uk/sasquatch/cgi-bin/foot.cgi) was run using default Workflow 3 settings (7-mer, propensity-based [Erythroid], exhaustive) on the web interface.

Analysis scripts for element accessibility are in *enhancer_accessibility*


**Allelic bias analysis**

[WASP](https://github.com/bmvdgeijn/WASP) used an adapted Snakemake pipeline. More information can be found in the CRAN package [documentation](https://cran.r-project.org/web/packages/coloc/index.html). Adapted scripts are in: *wasp_pipe_for_covid_gwas*


**3C analysis**

Domain inferral for probe design used the [3D genome browswer](http://3dgenome.fsm.northwestern.edu/index.html) and probes were designed using [Capsequm2](https://apps.molbiol.ox.ac.uk/CaptureC/cgi-bin/CapSequm.cgi).

NuTi/NG Capture-C replicate data were first processed with [CCSeqBasic5](https://github.com/Hughes-Genome-Group/CCseqBasicS) before being combined with [CaptureCompare](https://github.com/Hughes-Genome-Group/CaptureCompare). Run scripts and custom CaptureCompare scripts are in *capture_c_scripts*


MCC analysis codes are available for [academic use](https://process.innovation.ox.ac.uk/software/p/16529a/micro-capture-c-academic/1) with input files in *mcc_scripts*.

[LanceOtron](https://github.com/Hughes-Genome-Group/Lanceotron-User-Docs)

**Expression analysis**

GTEx [Calculator](https://www.gtexportal.org/home/testyourown)


**Colocalisation analysis**

Colocalisation of eQTL and GWAS signals used [Coloc](https://github.com/chr1swallace/coloc). The runshell used is *coloc_analysis.R*


**Editing ICE analysis**

Editing efficiency calculation used the Synthego [ICE](https://ice.synthego.com/#/) webtool.


**ChIP/ATAC-seq mapping**

Reads were mapped using [NGseqBasic](https://github.com/Hughes-Genome-Group/NGseqBasic)


**Spatial Transcriptomics**

Modules and cell types were identified using [WGCNA](https://rdrr.io/cran/WGCNA/) and [Deconvolution](https://rdrr.io/bioc/SpatialDecon/src/R/package.R). For further information see [Cross et al](https://www.biorxiv.org/content/10.1101/2021.06.21.449178v1). Further analysis was implemented in R with: *Spatial_analysis.R* 
