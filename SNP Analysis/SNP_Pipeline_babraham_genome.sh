#!/bin/bash
#Initial STAR Pass

mkdir -p STAR_babraham

for i in `find . -name "*_R1_trim.fq.gz" | sed 's/_R1_trim.fq.gz//'`
	do
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i 1st Pass Alignment begun <<<<<<<<<<<<<<<<<<<<"
		echo ""
		STAR --genomeDir /mnt/Scratch/Swine_flu/Sus_scrofa/Babraham_STAR \
    --runMode alignReads \
		--readFilesIn $i\_R1_trim.fq.gz $i\_R2_trim.fq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix STAR_babraham/`basename $i` \
    --outSJfilterReads Unique \
		--outSAMtype BAM SortedByCoordinate --runThreadN 60
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i 1st Pass Alignment complete <<<<<<<<<<<<<<<<<<<<"
		echo ""
	done

#command; cat STAR_1st/*final* | grep -i Unique | mail -s "1st Pass Alignments Complete" a.paterson@bristol.ac.uk

mkdir -p STAR_Babraham_2nd

for i in `find . -name "*_R1_trim.fastq.gz" | sed 's/_R1_trim.fastq.gz//'`
	do
		base="$(basename "$i")"
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i 2nd Pass Alignment begun <<<<<<<<<<<<<<<<<<<<"
		echo ""
		STAR --runThreadN 100 \
		--genomeDir /mnt/Scratch/Swine_flu/Sus_scrofa/Babraham_STAR \
		--outFileNamePrefix STAR_Babraham_2nd/`basename $i` \
		--readFilesIn $i\_R1_trim.fastq.gz $i\_R2_trim.fastq.gz \
		--readFilesCommand zcat \
		--outSJfilterReads Unique \
		--outSAMstrandField intronMotif \
    --outSAMunmapped Within \
		--outFilterMultimapNmax 1 \
		--outSAMattributes NH HI AS NM MD \
		--outSAMattrRGline ID:$base SM:$base LB:$base\-RNAseq PU:$base\-202003 PL:illumina \
		-sjdbFileChrStartEnd STAR_babraham/*.out.tab \
		--outSAMtype BAM SortedByCoordinate
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i 2nd Pass Alignment Complete <<<<<<<<<<<<<<<<<<<<"
		echo ""
	done

# Merge Alignments - Picard
# Unneccessary - STAR reports unaligned within same bam file output.

# Mark Duplicates - MarkDuplicatesSpark

for i in `find STAR_2nd_v2/. -name "*bam" | sed 's/.bam//'`
  do
		base="$(basename "$i")"
		gatk MarkDuplicatesSpark \
    -I $i\.bam \
    -O STAR_2nd_v2/$base\_marked_duplicates.bam \
    -M STAR_2nd_v2/$base\_marked_dup_metrics.txt
  done

#Dict file AND index required
# samtools faidx /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna
# java -jar picard.jar CreateSequenceDictionary R=/mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna O=/mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.dict

# Reformat RNAseq alignments that span introns for HaplotypeCaller.
# This step splits reads with N in the cigar into multiple supplementary alignments and hard clips mismatching overhangs.
# By default this step also reassigns mapping qualities for good alignments to match DNA conventions

for i in `find STAR_2nd_v2/. -name "*_marked_duplicates.bam" | sed 's/_marked_duplicates.bam//'`
  do
		base="$(basename "$i")"
    gatk SplitNCigarReads \
    -R /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna \
    -I STAR_2nd_v2/$base\_marked_duplicates.bam \
    -O STAR_2nd_v2/$base\_splitNcigarreads.bam
  done

# Variant Calling
for i in `find STAR_2nd_v2 -name "*_splitNcigarreads.bam" | sed 's/_splitNcigarreads.bam//'`
	do
		base="$(basename "$i")"
		echo $base
		echo STAR_2nd_v2/$base\_splitNcigarreads.bam
		gatk --java-options "-Xmx500g" HaplotypeCaller \
		-ERC GVCF \
  	-R /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna \
  	-I STAR_2nd_v2/$base\_splitNcigarreads.bam \
  	-O STAR_2nd_v2/$base\_output.g.vcf.gz
	done

# Variant Calling - 1st half
	for i in `find STAR_2nd -name "[A-D]*_splitNcigarreads.bam" | sed 's/_splitNcigarreads.bam//'`
		do
			base="$(basename "$i")"
			echo $base
			echo STAR_2nd_v2/$base\_splitNcigarreads.bam
			gatk --java-options "-Xmx400g" HaplotypeCaller \
			-ERC GVCF \
	  	-R /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna \
	  	-I STAR_2nd/$base\_splitNcigarreads.bam \
	  	-O STAR_2nd/$base\_output.g.vcf.gz
		done
# Variant Calling - 2nd half

		for i in `find STAR_2nd -name "[E-H]*_splitNcigarreads.bam" | sed 's/_splitNcigarreads.bam//'`
			do
				base="$(basename "$i")"
				echo $base
				echo STAR_2nd_v2/$base\_splitNcigarreads.bam
				gatk --java-options "-Xmx400g" HaplotypeCaller \
				-ERC GVCF \
				-R /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna \
				-I STAR_2nd/$base\_splitNcigarreads.bam \
				-O STAR_2nd/$base\_output.g.vcf.gz
			done


# Create intervals.list required for GenomicsDBImport
grep "^>" /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna | sed 's/ .*//g' | sed 's/>//g' > intervals.list


# Sample merging
#Babrahams
gatk GenomicsDBImport \
-V H05-B06-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V D06-B08-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V D08-B106-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V B06-B160-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V D05-B166-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V A08-G04-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V C08-G12-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V A06-G22-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V F07-G24-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V F08-G28-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V G05-G29-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V G07-G30-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V B08-G34-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V B05-G37-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V B07-G40-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V D07-G47-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V C06-P013-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V F06-P107-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V E06-P120-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V E05-P156-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V C05-P162-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V H04-P164-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V F05-P236-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V C07-R06-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V H07-R111-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V A05-R71-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V G06-R73-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V H06-R92-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V E07-R97-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V E08-R98-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V A07-Y204-ATAligned.sortedByCoord.out_output.g.vcf.gz \
--genomicsdb-workspace-path babraham_snps \
--tmp-dir=babraham_snps_temp \
-L intervals.list

#Commercial
gatk GenomicsDBImport \
-V E01-B22-ATAligned.sortedByCoord.out_output.g.vcf.gz \
-V H02-B36-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V A01-B42-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V G01-B72-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V C02-G02-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V G02-G104-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V A02-G105-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V D01-G107-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V G04-G114-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V E03-G115-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V A04-G116-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V C04-G117-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V D04-P150-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V H03-P224-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V C03-P349-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V H01-R09-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V D03-R146-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V E04-R147-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V F04-R148-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V F03-R149-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V C01-R49-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V E02-R59-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V B01-R69-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V F02-Y201-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V B04-Y272-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V G03-Y283-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V B02-Y355-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V F01-Y361-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V B03-Y385-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V A03-Y386-ATAligned.sortedByCoord.out_output.g.vcf.gz\
-V D02-Y397-ATAligned.sortedByCoord.out_output.g.vcf.gz\
--genomicsdb-workspace-path commercial_snps \
-L intervals.list \
--tmp-dir=commercial_snps_temp


# Joint-Call Cohort
gatk --java-options "-Xmx400g" GenotypeGVCFs \
	-R /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna \
	-V gendb://babraham_snps \
	-O babraham.vcf.gz

gatk --java-options "-Xmx400g" GenotypeGVCFs \
	-R /mnt/Scratch/Swine_flu/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna \
	-V gendb://commercial_snps \
	-O commercial.vcf.gz


gunzip *.gz

# Output SNPs only
gatk SelectVariants \
    -V babraham.vcf \
    -select-type SNP \
    -O babraham_snps.vcf

gatk SelectVariants \
  	-V commercial.vcf \
		-select-type SNP \
		-O commercial_snps.vcf


# Filtering - as per https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
gatk VariantFiltration \
   -V babraham_snps.vcf \
	 -filter "QD < 2.0" --filter-name "QD2" \
	 -filter "QUAL < 30.0" --filter-name "QUAL30" \
	 -filter "SOR > 3.0" --filter-name "SOR3" \
	 -filter "FS > 60.0" --filter-name "FS60" \
	 -filter "MQ < 40.0" --filter-name "MQ40" \
	 -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	 -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	 -O babraham_snps_filtered.vcf

gatk VariantFiltration \
  -V commercial_snps.vcf \
 	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "SOR > 3.0" --filter-name "SOR3" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	-O commercial_snps_filtered.vcf

# Convert NCBI chromosome accessions to chromosome numbers

sed	's/NC_010443\.5/1/g
	s/NC_010444\.4/2/g
	s/NC_010445\.4/3/g
	s/NC_010446\.5/4/g
	s/NC_010447\.5/5/g
	s/NC_010448\.4/6/g
	s/NC_010449\.5/7/g
	s/NC_010450\.4/8/g
	s/NC_010451\.4/9/g
	s/NC_010452\.4/10/g
	s/NC_010453\.5/11/g
	s/NC_010454\.4/12/g
	s/NC_010455\.5/13/g
	s/NC_010456\.5/14/g
	s/NC_010457\.5/15/g
	s/NC_010458\.4/16/g
	s/NC_010459\.5/17/g
	s/NC_010460\.4/18/g
	s/NC_010461\.5/X/g
	s/NC_010462\.3/Y/g' \
	$1 > $2


# Sample inspecting

# Total number of SNPs
bcftools view -v snps commercial_snps_filtered_chrnum.vcf.gz | grep -v "^#" | wc -l
860273

bcftools view -v snps babraham_snps_filtered_chrnum.vcf.gz | grep -v "^#" | wc -l
619099

# Total number of unique positions
bcftools view -v commercial_snps_filtered_chrnum.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
857427

bcftools view -v snps babraham_snps_filtered_chrnum.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
617633
