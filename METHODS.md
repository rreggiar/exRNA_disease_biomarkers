
Building new *salmon* references for `v1.6.0` and gencode `v38`:
```
./00c_linkProjectDir.sh /public/groups/kimlab/seqData/2020-09-30_exoRNABiomarkersPancAndCovid/bioIvt_panc
exRNA_disease_biomarkers
/public/groups/kimlab/seqData/2020-09-30_exoRNABiomarkersPancAndCovid/bioIvt_panc/ROSTER.csv
/public/groups/kimlab/exRNA_disease_biomarkers/data/input_data/bioIvt_panc
```

Linking in raw data:
```nohup ./00a_buildGencodeDirectory.sh 38 2>&1 > ../tmp/logs/00a_salmon_1.6.0_gencode_38_log.txt &```
