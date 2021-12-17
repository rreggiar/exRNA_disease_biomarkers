#!/bin/bash

INPUTDIR='/public/groups/kimlab/exoRNA-biomarkers-panc/data/'
TETXINDEX='/public/groups/kimlab/indexes/sel.aln.gen.34.ucsc.rmsk.index.salmon.1.2.1/'
GENINDEX='/public/groups/kimlab/indexes/sel.aln.gen.34.index.salmon.1.2.1/'
HG38='/public/groups/kimlab/genomes.annotations/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
GENOMEDIR='/public/groups/kimlab/genomes.annotations/HG.38.w.te.cons/'
ADAPTERS='/public/groups/kimlab/genomes.annotations/adapters'
TEOUTDIR=te.locus.gencode.34.salmon.out
GENOUTDIR=gencode.34.salmon.out

STAR --version
salmon -v
trimmomatic -version

cd $INPUTDIR

for SET in $PWD/*; do

  if [ $(basename $SET) = 'panc.new.set' ]; then

    LIB=TruSeq3-PE.fa

  else

    LIB=NexteraPE-PE.fa

  fi

  for INDEX in $TETXINDEX $GENINDEX; do

    cd $SET

    for SAMPLE in $PWD/*; do
        
        echo $SAMPLE

	      cd $SAMPLE

	      NAME=$(basename $SAMPLE)

	      if [ ! -f *_paired.fq.gz ]; then
	
		      read1=*_R1_001.fastq.gz
		      read2=*_R2_001.fastq.gz
	        echo "Read 1:" ${read1}
		      echo "Read 2:" ${read2}
          echo "fastq files must be trimmed"
          echo 'Trimming Now..'
          trimmomatic PE ${read1} ${read2} output_forward_paired.fq.gz \
            output_forward_unpaired.fq.gz output_reverse_paired.fq.gz \
            output_reverse_unpaired.fq.gz \
            ILLUMINACLIP:$ADAPTERS/$LIB:1:30:10:4:true \
            LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 MINLEN:36
        fi

        echo "setting index"

        if [ $INDEX = $TETXINDEX ]; then
        
          OUTDIR=$TEOUTDIR

        else

          OUTDIR=$GENOUTDIR

        fi

	      if [ -f $SAMPLE/$OUTDIR/quant.sf ]; then
		      echo "SALMON HAS ALREADY BEEN RUN"

	      else
		      echo "SALMON HAS NOT BEEN RUN, NOW RUNNING" ${read1}
      
		
		      mkdir $SAMPLE/$OUTDIR/
		
          trim1=$SAMPLE/*output_forward_paired.fq.gz
          trim2=$SAMPLE/*output_reverse_paired.fq.gz

		      salmon quant -i ${INDEX} \
			      --libType A \
			      -1 ${trim1}\
			      -2 ${trim2}\
			      -p 32\
			      --validateMappings \
			      --gcBias \
			      --seqBias \
            --numBootstraps 10 \
            --recoverOrphans \
            --rangeFactorizationBins 4 \
			      --output $SAMPLE/$OUTDIR/

        fi
    
        echo "SALMON COMPLETE, CHECKING STAR"

        trim1=$SAMPLE/output_forward_paired.fq.gz
        trim2=$SAMPLE/output_reverse_paired.fq.gz

        if [ ! -f $SAMPLE/star.out/Aligned.out.sam ]; then

          echo "RUNNING 2 PASS STAR ON" $SAMPLE
          firstRunDir=$SAMPLE/star.out/
          mkdir $firstRunDir
          cd $firstRunDir

          STAR --genomeDir $GENOMEDIR \
            --readFilesIn <(gunzip -c $trim1) <(gunzip -c $trim2) \
            --runThreadN 32 \
            --outMultimapperOrder Random \
            --outFilterMultimapNmax 50 \
            --outReadsUnmapped Fastx
        else
          
          echo "STAR IS COMPLETE"

        fi

    done

  echo "DONE WITH FIRST PASS STAR AND SALMON"

  done

  cd $SET

  for SAMPLE in $PWD/*; do

    if [ ! -d $SAMPLE/star.out/pass.2 ]; then

      secondRunDir=$SAMPLE/star.out/pass.2/
      mkdir $secondRunDir


      trim1=$SAMPLE/output_forward_paired.fq.gz
      trim2=$SAMPLE/output_reverse_paired.fq.gz

      cd $secondRunDir

	    STAR --genomeDir $GENOMEDIR \
        --readFilesIn <(gunzip -c $trim1) <(gunzip -c $trim2) \
        --runThreadN 32 \
        --sjdbFileChrStartEnd $SET/*/star.out/SJ.out.tab \
        --outFilterMultimapNmax 50 \
        --outReadsUnmapped Fastx \
        --outMultimapperOrder Random \
        --outSAMtype BAM SortedByCoordinate
    fi

  done

done
