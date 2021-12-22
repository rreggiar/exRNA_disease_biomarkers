
Building new *salmon* references for `v1.6.0` and gencode `v38`:  
```
./00c_linkProjectDir.sh /public/groups/kimlab/seqData/2020-09-30_exoRNABiomarkersPancAndCovid/bioIvt_panc
exRNA_disease_biomarkers
/public/groups/kimlab/seqData/2020-09-30_exoRNABiomarkersPancAndCovid/bioIvt_panc/ROSTER.csv
/public/groups/kimlab/exRNA_disease_biomarkers/data/input_data/bioIvt_panc
```

Linking in raw data:  
```nohup ./00a_buildGencodeDirectory.sh 38 2>&1 > ../tmp/logs/00a_salmon_1.6.0_gencode_38_log.txt &```

Running *salmon* for gencode: [code](bin/01_wrapper.sh)  
```
nohup ./01_trimAndSalmon.sh /public/groups/kimlab/exRNA_disease_biomarkers/data/input_data/bioIvt_ctrl/ \
	exo \
	/public/groups/kimlab/indexes/sel.align.gencode.v38.salmon.v1.6.0.sidx \
	/public/groups/kimlab/exRNA_disease_biomarkers/data/output_data/rna_qc \
	2>&1 > /public/groups/kimlab/exRNA_disease_biomarkers/tmp/logs/bioIvt_ctrl_quant_log.txt &

....

```  

Generate a new *STAR* genome for `v2.7.9a`:  
```
nohup STAR --runThreadN 24 --runMode genomeGenerate \
	--genomeDir /public/groups/kimlab/genomes.annotations/hg38_star_2.7.9a \
	--genomeFastaFiles /public/groups/kimlab/genomes.annotations/GRCh38.p13.genome.fa \
	2>&1 > tmp/logs/STAR_genomeGenerate_2.7.9a_log.txt &
```  

Running *STAR*: [code](bin/02_wrapper.sh)  

Now also have *STAR* qc : [code](bin/02_wrapper_qc.sh)  
```
./02_wrapper_qc.sh '/public/groups/kimlab/exRNA_disease_biomarkers/'
```

