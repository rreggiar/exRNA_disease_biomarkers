#!/bin/bash

# array of jobs to execute the entire project's rna quant and qc workload

nohup ./01_trimAndSalmon.sh /public/groups/kimlab/exRNA_disease_biomarkers/data/input_data/bioIvt_ctrl/ \
	exo \
	/public/groups/kimlab/indexes/sel.align.gencode.v38.salmon.v1.6.0.sidx \
	/public/groups/kimlab/exRNA_disease_biomarkers/data/output_data/rna_qc \
	2>&1 > /public/groups/kimlab/exRNA_disease_biomarkers/tmp/logs/bioIvt_ctrl_quant_log.txt &

nohup ./01_trimAndSalmon.sh /public/groups/kimlab/exRNA_disease_biomarkers/data/input_data/bioIvt_covid \
	exo \
	/public/groups/kimlab/indexes/sel.align.gencode.v38.salmon.v1.6.0.sidx \
	/public/groups/kimlab/exRNA_disease_biomarkers/data/output_data/rna_qc \
	2>&1 > /public/groups/kimlab/exRNA_disease_biomarkers/tmp/logs/bioIvt_covid_quant_log.txt &

nohup ./01_trimAndSalmon.sh /public/groups/kimlab/exRNA_disease_biomarkers/data/input_data/bioIvt_luad \
	exo \
	/public/groups/kimlab/indexes/sel.align.gencode.v38.salmon.v1.6.0.sidx \
	/public/groups/kimlab/exRNA_disease_biomarkers/data/output_data/rna_qc \
	2>&1 > /public/groups/kimlab/exRNA_disease_biomarkers/tmp/logs/bioIvt_luad_quant_log.txt &


nohup ./01_trimAndSalmon.sh /public/groups/kimlab/exRNA_disease_biomarkers/data/input_data/bioIvt_panc \
	exo \
	/public/groups/kimlab/indexes/sel.align.gencode.v38.salmon.v1.6.0.sidx \
	/public/groups/kimlab/exRNA_disease_biomarkers/data/output_data/rna_qc \
	2>&1 > /public/groups/kimlab/exRNA_disease_biomarkers/tmp/logs/bioIvt_panc_quant_log.txt &


nohup ./01_trimAndSalmon.sh /public/groups/kimlab/exRNA_disease_biomarkers/data/input_data/panc1_intra/ \
	intra \
	/public/groups/kimlab/indexes/sel.align.gencode.v38.salmon.v1.6.0.sidx \
	/public/groups/kimlab/exRNA_disease_biomarkers/data/output_data/rna_qc \
	2>&1 > /public/groups/kimlab/exRNA_disease_biomarkers/tmp/logs/panc1_intra_quant_log.txt &


nohup ./01_trimAndSalmon.sh /public/groups/kimlab/exRNA_disease_biomarkers/data/input_data/pittsburgh_ctrl \
	exo \
	/public/groups/kimlab/indexes/sel.align.gencode.v38.salmon.v1.6.0.sidx \
	/public/groups/kimlab/exRNA_disease_biomarkers/data/output_data/rna_qc \
	2>&1 > /public/groups/kimlab/exRNA_disease_biomarkers/tmp/logs/pittsburgh_ctrl_quant_log.txt &

