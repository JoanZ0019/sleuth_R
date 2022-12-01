########### kallisto command line: ###########
#####Bulid an index of transcript.fats file # CODE： kallisto index -i anjou_transcripts.idx anjou_transcripts.fasta
##### Quantify expression #CODE：kallisto quant -t 4 -i anjou_transcripts.idx -o output ./RNA_seq/JLHW_RNAseq_0021_Expanding_Leaves_GATTCTGC_Pear_dAnjou_I1126_L1_R1.fastq.gz ./RNA_seq/JLHW_RNAseq_0021_Expanding_Leaves_GATTCTGC_Pear_dAnjou_I1126_L1_R2.fastq.gz
##### Bootstrap #CODE: kallisto quant -i ./anjou_transcripts.idx -t 4 -b 10 ./RNA_seq/JLHU_RNAseq_0016_Budding_Leaves_AGCCTCAT_Pear_dAnjou_I1126_L1_R1.fastq.gz ./RNA_seq/JLHU_RNAseq_0016_Budding_Leaves_AGCCTCAT_Pear_dAnjou_I1126_L1_R2.fastq.gz -o budding_leaves_bootstrap_output 
#___________________________________________________________________________________________________________________________________________________________________________

setwd("C:/Users/qzz0019/Box/Courses in Auburn University/2022 Fall/kallisto_OUTPUT")
## Step 1: install "sleuth" 
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()
BiocManager::install("devtools")    # only if devtools not yet installed
BiocManager::install("pachterlab/sleuth")
library("devtools")
library("sleuth")
## Load sleuth and accessory libraries
require("sleuth")
packageVersion("sleuth") #‘0.30.1’
# install.packages("gridExtra")
# install.packages("cowplot")
library("gridExtra")
library("cowplot")
# BiocManager::install("biomaRt") #to allow us to pull in recognizable gene names from a database.
library("biomaRt")
#---------------------------------------------------------------------------------------------------

## ________________Step 2: Load experimental design and label kallisto output with metabdata________________

## Locate sample names and describe the experimental design
#provide Sleuth with sample names
sample_id <- dir(file.path("C:/Users/qzz0019/Box/Courses in Auburn University/2022 Fall/kallisto_OUTPUT/kallisto_result"))
sample_id
#get the file paths to results files
kal_dirs<- file.path("C:/Users/qzz0019/Box/Courses in Auburn University/2022 Fall/kallisto_OUTPUT/kallisto_result", sample_id)
kal_dirs
#need a table that provides more meaningful names for describing the experiment, group samples as "replicates",like fruit, leaves, buds.
s2c <- read.table(file.path("C:/Users/qzz0019/Box/Courses in Auburn University/2022 Fall/kallisto_OUTPUT/experimental_design.txt"),  header = TRUE,   stringsAsFactors = FALSE)
s2c
#add file path to the table 
s2c <- dplyr::mutate(s2c,path = kal_dirs)
s2c
##read in the transcript-gene mapping file (extract from anjou_pros.fasta, the first column as transcript ids (the fasta headers in the protein fasta file without the > character) and the second column as the gene id (the transcript id without the period number at the end)
t2g <- read.table(file.path("C:/Users/qzz0019/Box/Courses in Auburn University/2022 Fall/kallisto_OUTPUT/t2g.txt"), header = TRUE, stringsAsFactors=FALSE)
#t2g
# ----------------------------------------------------------------------------------------------------------

## ______________Step 3: Construct the "Sleuth Objoect_____________" 

so <- sleuth_prep(s2c, ~tissue, target_mapping = t2g, extra_bootstrap_summary = TRUE)
#To generate bootstrap 
#CODE: kallisto quant -i ./anjou_transcripts.idx -t 4 -b 10 ./RNA_seq/JLHU_RNAseq_0016_Budding_Leaves_AGCCTCAT_Pear_dAnjou_I1126_L1_R1.fastq.gz ./RNA_seq/JLHU_RNAseq_0016_Budding_Leaves_AGCCTCAT_Pear_dAnjou_I1126_L1_R2.fastq.gz -o budding_leaves_bootstrap_output 

##  Fit the linear model and test for one of the model coefficients
 #1. Fit the full model
so <- sleuth_fit(so)
 #2. Fit to a "reduced" model that presumes abundances are equal in the two conditions
so <- sleuth_fit(so,~1, 'reduced')
models(so)
 #3. Sleuth will perform LRT to identify transcripts with a significantly better fit with the “full” model
so = sleuth_lrt(so, null_model = "reduced", alt_model = "full")

#Retrieve results of the test
res = sleuth_results(so, test = "reduced:full", test_type = "lrt", show_all = TRUE)
head(res)

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20) #The table shown above displays the top 20 significant genes with a (Benjamini-Hochberg multiple testing corrected) q-value <= 0.05.
#-------------------------------------------------------------------------------------------------

## ______________Step 4: Visualize results_____________" 

## <<<<<<<<<<   Plot clustered transcript HEATMAP    >>>>>>>>>>>>>>>
# Create expression HEATMAP for top 40 genes with lowest q-values
png("DGE_sleuth_transcript_heatmap.png", width=7, height=7, units = "in", res = 300)
plot_transcript_heatmap(so, cluster_transcripts = TRUE,
                        transcripts = subset(res[order(res$qval),], qval < 0.05)$target_id[1:40]) # top40
dev.off()

## <<<<<<<<<<   PCA" Principal component plot of the samples   >>>>>>>>>>
png("DGE_sleuth_PCA.png", width=7, height=7, units = "in", res = 300)
plot_pca(so, color_by = "tissue", text_labels = TRUE, units = "tpm")
dev.off()


##The easiest way to view and interact with the results is to generate the sleuth live site that allows for exploratory data analysis
sleuth_live(so)





