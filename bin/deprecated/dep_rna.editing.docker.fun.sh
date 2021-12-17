
for BEDFILE in sleuth.alu.bed sleuth.line1.bed sleuth.mer.bed sleuth.ltr.bed; do
#for BEDFILE in gencode.as1.bed kras.alu.bed kras.mer.bed sz.sx.sx1.y.bed; do
  docker run \
    -u $(id -u ${USER}):$(id -g ${USER})\
    -v '/public/groups/kimlab/aale.kras/data/single.cell.rna.seq/output/cluster.bam.files/':'/data/input_files'\
    -v '/public/groups/kimlab/aale.kras/data/single.cell.rna.seq/output/subset.rna.editing.out/'$BEDFILE:'/data/output_dir'\
    rna_editing:1.0 RNAEditingIndex -d '/data/input_files/'\
    -l '/data/output_dir/log.files' -o '/data/output_dir/cmpileups'\
    -os '/data/output_dir/summary_dir'\
    -f '.bam'\
    --genome hg38 --verbose -rb '/data/input_files/'$BEDFILE
done

    #-v '/public/groups/kimlab/aale.kras/data/single.cell.rna.seq/output/rna.editing.out/hg38.out':'/data/output_dir' \

