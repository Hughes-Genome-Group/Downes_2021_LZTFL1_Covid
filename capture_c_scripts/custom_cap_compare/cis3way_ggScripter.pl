#!/usr/bin/perl
use strict;
use Cwd;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;

# This script is for bining and windowing capture C data. Requires normalised unionised files for a single viewpoint from
# 2 different samples/tissues. Single imput parameter file containing exclusion regions, visualisation region, and bin/window sizes.
# Output is a file of parameters for each viewpoint, and tab delimited multibedgraph of bins and windows

# Script adapted from Lar's manual scripts to be applicapble to multiple loci and more user friendly
# R scripts developed by Ron Schwessinger
### (C) Damien Downes 16th May 2018.


&GetOptions
(
    "viewpoints=s"=>\ my $viewPoints,     # --viewpoints      VIEWPOINT	CHR VP_START VP_STOP EXCLSTART EXCLSTOP REGIONSTART REGIONSTOP BINSIZE WINDOWSIZE
	"annotation=s"=> \ my $annotation,    # --annotation 	  Full path to Bed file with genes or other features.
);

my $current_directory = cwd;


open (VIEWPOINTS, $viewPoints) or die "Cannot open viewpoint list\n";
while (my $viewpoint = <VIEWPOINTS>)
	{
		chomp($viewpoint);
		my ($viewID, $view_chr, $vp_start, $vp_stop, $excl_start, $excl_stop, $region_start, $region_stop, $bin_size, $window_size) = split(' ', $viewpoint);
        
        open (OUTFH, ">$viewID\_plotting.R") or die "Cannot open file $viewID\_plotting.R\n\n";
		

		print OUTFH "# Set Working Directory and Load Libraries=====================================\n";		
        print OUTFH "setwd(\"$current_directory\")\n\n";
		
		print OUTFH "library(tidyverse)\n";
		print OUTFH "library(cowplot)\n";
		print OUTFH "library(RColorBrewer)\n\n";
		
		print OUTFH "# Plotting Theme ==============================================================\n";
		print OUTFH "science_theme <- theme(\n";
		print OUTFH "  panel.grid.major = element_line(size = 0.5, color = \"grey\"),\n";
		print OUTFH "  panel.grid.minor = element_blank(),\n";
		print OUTFH "  plot.title = element_text(hjust = 0.5),\n";
		print OUTFH "  text = element_text(size = 14),";
		print OUTFH "  axis.line = element_line(color=\"black\", size = 0.7),\n";
		print OUTFH "  axis.line.x = element_line(color=\"black\", size = 0.7),\n";
		print OUTFH "  axis.line.y = element_line(color=\"black\", size = 0.7),\n";
		print OUTFH "  # plot.margin = unit(c(0.7,0.7,0.7,0.7), \"lines\"),\n";
		print OUTFH "  panel.border=element_blank(),\n";
		print OUTFH "  strip.background = element_blank()\n";
		print OUTFH ")\n\n";
		
		print OUTFH "# Helper Function =============================================================\n";
		print OUTFH "IntersectBedDataframe <- function(df, chr, start, end){\n";
		print OUTFH "  # Wrapper to subset/intersect a given BED format df with a region of interest\n";
		print OUTFH "  #\n";
		print OUTFH "  # Input:\n";
		print OUTFH "  #   df: bedlike 4 or more column df\n";
		print OUTFH "  #   chr: chromsome of interest\n";
		print OUTFH "  #   start: start coord of interest\n";
		print OUTFH "  #   end: end coord of interest\n";
  
		print OUTFH "  df <- df[df[,1] == chr,] #chr\n";
  
		print OUTFH "  df <- df[\n";
		print OUTFH "    (df[, 2] >= start & df[, 2] < end) |\n";
		print OUTFH "      (df[, 3] > start & df[, 3] <= end) |\n";
		print OUTFH "      (df[, 2] <= start & df[, 3] >= end),] # coords\n";
  
		print OUTFH "  return(df)\n";
  
		print OUTFH "}		\n\n";
		
		
		print OUTFH "# Load Data and Parameter File ================================================\n";
		
		print OUTFH "data.file <- \"$viewID\_window.tab\"\n";
		print OUTFH "parameters.file <- \"$viewPoints\"\n";
		print OUTFH "id <- \"$viewID\"\n\n";	
		
		print OUTFH "data <- as_tibble(read.table(data.file, header=T))\n";
		print OUTFH "parameters <- as_tibble(read.table(parameters.file))\n";
		print OUTFH "names(parameters) <- c(\"Viewpoint\", \"Chr\", \"Frag_start\", \"Frag_stop\", \"Exclusion_Start\", \"Exclusion_Stop\", \"Plot_Region_Start\", \"Plot_Region_Stop\", \"Bin_Size\", \"Window_size\")\n\n";
		

		print OUTFH "# Modify and select from parameters file =======================================\n";
		print OUTFH "# select viewpoint from parameters file\n";
		print OUTFH "viewp <- parameters %>% filter(Viewpoint == id)\n";
		
		print OUTFH "# select exclusion fragments within plotting region\n";
		print OUTFH "parameters <- parameters %>%\n";
		print OUTFH "  filter(Frag_stop >= viewp\$Plot_Region_Start & Frag_start <= viewp\$Plot_Region_Stop)\n";
		print OUTFH "# add window size to exclusion zone only for viewpoint\n";
		print OUTFH "parameters <- parameters %>%\n";
		print OUTFH "  mutate(Exclusion_Start = if_else(Viewpoint == id, (Exclusion_Start - Window_size), Exclusion_Start)) %>%\n";
		print OUTFH "  mutate(Exclusion_Stop = if_else(Viewpoint == id, (Exclusion_Stop + Window_size), Exclusion_Stop))\n\n";
		
		print OUTFH "to.exclude <- parameters %>%\n";
		print OUTFH "	filter(Viewpoint == id)\n\n";

			#		- Could hard code in HbaCombined and HbbCombined for exclusion of both windows.
		
		print OUTFH "# Load Data for Gene intersect =================================================\n";
		print OUTFH "genes.file <- \"$annotation\"\n";
		print OUTFH "known.genes <- read.table(genes.file)\n";

		print OUTFH "# Intersect Genes with plot frame ==============================================\n";
		print OUTFH "genes <- as_tibble(IntersectBedDataframe(known.genes, chr=as.character(viewp\$Chr), start=viewp\$Plot_Region_Start, end=viewp\$Plot_Region_Stop))\n";
		print OUTFH "names(genes) <- c(\"chr\", \"start\", \"end\", \"name\")\n";
		print OUTFH "# some formating for plotting\n";
		print OUTFH "genes <- genes %>% \n";
		print OUTFH "  mutate(annotation = \"Genes\") %>%\n";
		print OUTFH "  mutate(x = start + (end - start)/2) %>%\n";
		print OUTFH "  mutate(y = 2)\n";
		print OUTFH "genes[seq(1,nrow(genes), 2),\"y\"] <- 0 \n\n";
		

		print OUTFH "# Gather Data =================================================================\n";
		print OUTFH "# 1 gather in long format\n";
		print OUTFH "d <- data %>%\n";
		print OUTFH "  gather(key, value, -c(BinNum, Chr, Start, Stop))\n";
		print OUTFH "# 2 extract condition and replicate\n";
		print OUTFH "d <- d %>%\n";
		print OUTFH "  mutate(condition = sub('_.+', \"\", key, perl=T)) %>%\n";
		print OUTFH "  mutate(replicate = sub('.+_', \"\", key, perl=T)) %>%\n";
		print OUTFH "  mutate(pos = Start + (Stop - Start)/2)\n\n";

		
		print OUTFH "# Calculate Mean and STDEV ====================================================\n";
		print OUTFH "d <- d %>%\n";
		print OUTFH "  group_by(BinNum, Chr, Start, Stop, pos, condition) %>%\n";
		print OUTFH " summarize(mean = mean(value), sd = sd(value)) %>%\n";
 		print OUTFH " ungroup()\n\n";
		

		print OUTFH "# Plot Means ± StDev ==========================================================\n";
		print OUTFH "p <- ggplot(d, aes(x=pos, y=mean, col = condition, fill = condition)) + \n";
		print OUTFH "  geom_ribbon(inherit.aes = F, aes(x=pos, ymin=mean-sd, ymax=mean+sd, fill=condition), alpha=.35) + \n";
		print OUTFH "  geom_line() + \n";
		print OUTFH "  scale_fill_brewer(palette = \"Set1\") +  \n";
		print OUTFH "  scale_colour_brewer(palette = \"Set1\") +  \n";
		print OUTFH "  geom_rect(data=to.exclude, inherit.aes=F, aes(xmin=Exclusion_Start, xmax=Exclusion_Stop), ymin=-25, ymax=max(d\$mean)+max(d\$sd), col=\"lightgrey\", fill=\"white\") + \n";
		print OUTFH "  geom_rect(data=parameters, inherit.aes=F, aes(xmin=Frag_start, xmax=Frag_stop), ymin=-25, ymax=max(d\$mean)+max(d\$sd), col=\"grey\", fill=\"grey\") + \n";
		print OUTFH "  labs(x=\"position\", y=\"mean interaction\") + \n";
		print OUTFH "  coord_cartesian(xlim=c(viewp\$Plot_Region_Start, viewp\$Plot_Region_Stop)) + \n";
		print OUTFH "  ggtitle(id) +\n"; 
		print OUTFH "  theme_bw() + science_theme +  \n";
		print OUTFH "  theme(panel.grid = element_blank(),\n";
		print OUTFH "      axis.line.x = element_blank(),\n";
		print OUTFH "      axis.text.x = element_blank(),\n";
		print OUTFH "      axis.ticks.x = element_blank(),\n";
		print OUTFH "      axis.title.x = element_blank()\n";
		print OUTFH "      )	\n\n";
		
		print OUTFH "# Plot of Genes ===============================================================\n";
		print OUTFH "g <- ggplot(genes, aes(col=annotation)) +\n";
		print OUTFH "  geom_hline(yintercept = 1, col = brewer.pal(3, \"Set1\")[2]) +\n";
		print OUTFH "  geom_rect(aes(xmin=start, xmax=end, ymin=0.55, ymax=1.55), fill=\"white\") +\n";
		print OUTFH "  scale_colour_manual(values = brewer.pal(3, \"Set1\")[2]) + \n";
		print OUTFH "  geom_text(aes(label = name, x = x, y = y), size=3.5, col=\"black\", fontface=2) +\n";
		print OUTFH "  coord_cartesian(xlim=c(viewp\$Plot_Region_Start, viewp\$Plot_Region_Stop), ylim = c(-0.5, 2.5)) +\n";
		print OUTFH "  labs(x=\"position\", y=\"genes\") + \n";
		print OUTFH "  science_theme + \n";
		print OUTFH "  theme(\n";
		print OUTFH "    legend.position=\"None\",\n";
		print OUTFH "    axis.text.y=element_blank(), \n";
		print OUTFH "    axis.line.y=element_blank(),\n";
		print OUTFH "    axis.ticks.y = element_blank(),\n";
		print OUTFH "    panel.grid.major = element_blank(),\n";
		print OUTFH "    panel.border = element_blank(),\n";
		print OUTFH "    plot.margin = unit(c(0,1,0,1), \"lines\"))\n\n";

		print OUTFH "c <- plot_grid(p, g, align = \"v\", axis = \"lrtb\", nrow=2, rel_heights = c(5,1))\n\n";
		
		
		print OUTFH "ggsave(c, filename = \"$viewID\_plot.pdf\", width=12, height=6)\n";
		

        close OUTFH;
    }
close VIEWPOINTS;
exit;
