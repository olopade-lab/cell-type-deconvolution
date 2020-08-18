#!/bin/bash

python tpmProcessorPipeline.py /Users/arvindkrishnan/Downloads/salmon_merged_gene_tpm.csv all -o allGenes
python tpmProcessorPipeline.py /Users/arvindkrishnan/Downloads/salmon_merged_gene_tpm.csv quantiseq -o quantiseq
python tpmProcessorPipeline.py /Users/arvindkrishnan/Downloads/salmon_merged_gene_tpm.csv mcp_counter -o mcp_counter
python tpmProcessorPipeline.py /Users/arvindkrishnan/Downloads/salmon_merged_gene_tpm.csv cibersort -o cibersort
python tpmProcessorPipeline.py /Users/arvindkrishnan/Downloads/salmon_merged_gene_tpm.csv cibersort_abs -o cibersort_abs
python tpmProcessorPipeline.py /Users/arvindkrishnan/Downloads/salmon_merged_gene_tpm.csv timer -o timer
python tpmProcessorPipeline.py /Users/arvindkrishnan/Downloads/salmon_merged_gene_tpm.csv xcell -o xcell
python tpmProcessorPipeline.py /Users/arvindkrishnan/Downloads/salmon_merged_gene_tpm.csv epic -o epic