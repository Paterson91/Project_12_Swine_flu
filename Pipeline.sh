#!/bin/bash



#INDEXING SUS SCROFA GENOME (NCBI)

STAR --runThreadN 50 \
--runMode genomeGenerate \
--genomeDir /mnt/Scratch/Swine_flu/Sus_scrofa/STAR/ \
--genomeFastaFiles /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna \
--sjdbGTFfile /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.gff \
--sjdbOverhang 74 \
--limitGenomeGenerateRAM=500000000000

if [ $? -ne 0 ]
then
    command; echo "Indexing Error" | mail -s "Indexing Error" a.paterson@bristol.ac.uk
else
    command; echo "Indexing Complete" | mail -s "Indexing Complete" a.paterson@bristol.ac.uk
fi

#RUNNING ALIGNMENTS

mkdir -p STAR

for i in `find . -name "*_R1_trim.fastq.gz" | sed 's/_R1_trim.fastq.gz//'`
	do
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i 1st Pass Alignment begun <<<<<<<<<<<<<<<<<<<<"
		echo ""
		STAR --genomeDir /mnt/Scratch/Swine_flu/Sus_scrofa/STAR \
    --runMode alignReads \
		--readFilesIn $i\_R1_trim.fastq.gz $i\_R2_trim.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix STAR/`basename $i` \
    --outSJfilterReads Unique \
		--outSAMtype BAM SortedByCoordinate --runThreadN 100
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i 1st Pass Alignment complete <<<<<<<<<<<<<<<<<<<<"
		echo ""
	done

#Parsing GFF3 to GTF for featureCounts

gffread GCF_000003025.6_Sscrofa11.1_genomic.gff -T -F -o GCF_000003025.6_Sscrofa11.1_genomic_2.gtf
perl -ne 'chomp; @a=split/\t/; %h=split(/ /,$a[8]); $a[8]=join(" ",("gene_name",$h{"gene_name"},"transcript_id",$h{"transcript_id"})); \
print join("\t",@a),"\n";' GCF_000003025.6_Sscrofa11.1_genomic_2.gtf > GCF_000003025.6_Sscrofa11.1_genomic_3.gtf


#COUNTING FEATURES

# run featureCounts on ALL
featureCounts \
-p -g gene_name -t exon \
-a /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic_3.gtf \
-o Counts_stranded_gtf.gene -F GTF -T 64 \
-s 0 \
*.bam

#Tidy Output - Remove header, negative controls, strand information
tail -n +2 Counts_stranded_gtf.gene | cut -f1,7-14,16-23,25-32,34-41,43-50,52-66,68-74 | sed 's/-ATAligned.sortedByCoord.out.bam//g'> Counts_stranded_gtf.gene_tidied

#Vim to remove BAM extensions...


#CANT SED REMOVE BAM EXTENSION?
#tail -n +2 Counts_stranded_gtf.gene | cut -f1,7- | sed 's/-.*[^ ]//'> Counts_stranded_gtf.gene_tidied



#Differential Expression Analysis - DESeq2
#NOTE; Edit config.R to change options

Rscript DESeq2.R config.R
