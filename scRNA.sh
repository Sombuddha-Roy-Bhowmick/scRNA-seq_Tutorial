fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $threads --html vdj_v1_hs_nsclc_5gex.fastp.html --json vdj_v1_hs_nsclc_5gex.fastp.json -i *R1*.fastq.gz  -o fastp_output/vdj_v1_hs_nsclc_5gex.trimmed.R1.fastq -I *R2*.fastq.gz  -O fastp_output/vdj_v1_hs_nsclc_5gex.trimmed.R2.fastq

cellranger count \
	--transcriptome refdata-gex-GRCh38-2020-A `# GRCh38 reference genome` #curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz \
	--id vdj_v1_hs_nsclc_5gex `# Sample ID` \
        --sample vdj_v1_hs_nsclc_5gex `# Sample name` \
        --expect-cells=1000 `# Expected number of cells` \
        --localcores=64 `# Number of CPU cores` \
        --localmem=250 `# Amount of memory (GB)` \
        --fastqs fastp_output/

Rscript scRNA.R  #Analysis of scRNA data 
