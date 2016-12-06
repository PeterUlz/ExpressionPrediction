# Expression Prediction
Analyze gene expression based on coverage distribution around TSS from cfDNA

Basic Usage:
expression_prediction.py [-h] -fq FASTQ_FILE -s NAME -g {m,f} [-o OUTDIR] [-k] [-t THREADS] -cna CNA_FILE [-tmp TEMP_DIR] [-step START_STEP]

Note:
This calls every step of the below list. Copy number alterations (CNAs) need to be normalized against. This is done by specifying a list of 
copy-number states in the format <chrom><start><stop><log2-ratio>. For plasma-seq analyses this is the *.segments file.

*Attention: Large hg19 reference index files for BWA are not in this commit*

Also needed for analysis:
. java 
. R (Package e1071)

-) Step1 create directory and create MD5 file of input
-) Step2 Trim fastq
-) Step3 alignment and conversion to BAM
-) Step4 remove PCR duplicates
-) Step5 analyze TSS profile in housekeeping vs. unexpressed genes and plot in R
-) Step6 extract coverage parameters for expression prediction
-) Step7 expression prediction

