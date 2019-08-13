#!/bin/bash
# Usage: bash nastybugs3.sh id.txt ./hgDir2 ./cardgene2 ./cardsnp2 16 ./outDir2 ./inDir 
# Required: id.txt ( List of SRA Accession Numbers - One Acession number per line)

if [ "$#" -ne 7 ]; then
	echo "Wrong number of arguments. Please enter bash nastybugs3.sh [List of SRA runs] [Host Genome Directory] [CARD Gene Database Directory] [CARD Variant Database Directory] [Cores] [Output Directory] [astq Directory]"
	exit 1
fi

sraid=$(readlink -f "$1") 
###JD:removed space after = in sraid assignment
hostGen=$(readlink -f "$2")
cardgene=$(readlink -f "$3")
cardsnp=$(readlink -f "$4")
cores="$5"
outdir=$(readlink -f "$6")
indir=$(readlink -f "$7")

## Getting host genome from NCBI and create BLAST db 

###JD:added if statement to check whether database has already been downloaded 
#if [ ! -e "$hostGen/GRCh37_latest_genomic.fna" ]; then
#  printf "Getting host genome and create BLAST databases\n"
#  date
#  mkdir $hostGen
#  wget -P $hostGen ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
#  gunzip -v $hostGen/GRCh37_latest_genomic.fna.gz
#  printf "Downloaded BLAST databases, start creating Bowtie index\n"
#  date
#  bowtie2-build $hostGen/GRCh37_latest_genomic.fna $hostGen/hg19 --threads 16
#  printf "Bowtie index created\n"
#  date
#fi

## Getting CARD Gene and SNP databases and create BLAST db
## https://card.mcmaster.ca/download
###JD:added if statement to check card db already exists
#if [ ! -e "$cardgene/cardgenedb.nhr" ]; then
#  printf "Getting CARD databases\n"
#  date
#  mkdir $cardgene $cardsnp
  ###JD:had to add --no-check-certificate in order to download - not sure if this is a good solution
#  wget --no-check-certificate -O $cardgene/broadstreet.tar.gz https://card.mcmaster.ca/latest/data
  ###JD:changed download location to latest rather than specific version
#  tar xvf $cardgene/broadstreet.tar.gz -C $cardgene/
#  printf "Downloaded CARD databases, start building CARD Bowtie indecies\n"
#  date
#  makeblastdb -in $cardgene/nucleotide_fasta_protein_homolog_model.fasta -parse_seqids -dbtype nucl -out $cardgene/cardgenedb
#  makeblastdb -in $cardgene/nucleotide_fasta_protein_homolog_model.fasta -dbtype nucl -out $cardgene/cardgenedb
#  bowtie2-build $cardgene/nucleotide_fasta_protein_homolog_model.fasta $cardgene/cardgenebt --threads 16
#  cp $cardgene/nucleotide_fasta_protein_variant_model.fasta $cardsnp/.
#  bowtie2-build $cardsnp/nucleotide_fasta_protein_variant_model.fasta $cardsnp/cardsnpbt --threads 16
#  printf "Done CARD Bowtie indecies\n"
#  date
#  makeblastdb -in $cardsnp/nucleotide_fasta_protein_variant_model.fasta -parse_seqids -dbtype nucl -out $cardsnp/cardsnpdb
#  makeblastdb -in $cardsnp/nucleotide_fasta_protein_variant_model.fasta -dbtype nucl -out $cardsnp/cardsnpdb
#fi

#mkdir $outdir

###JD:added prompt notification for user peace of mind
printf "\nNatsybugs3 processing reads...\n"

for sra in $(cat "$1"); do

#  printf "processing %s\n" "$sra"
#  date
#  fastq-dump --split-3 $sra -O $indir
  # First we align to a host so we can subtract host reads
###  magicblast -sra $sra -db $hostGen/hg19 -num_threads $cores -score 50 -penalty -3 -out $outdir/$sra.human.sam
#  printf "Dumped %s\n" "$sra"
#  date
#  

#  printf "aligned and unmapped output%s\n" "$sra"
#  date
  # Convert magicblast output sam file to FASTA
  # Getting unmapped reads (-f 4).  For mapped reads, use flag (-F 4)
  #if [[ (-s $outdir/$sra.human.sam) ]]; then
  #    samtools fasta -f 4 $outdir/$sra.human.sam -1 $outdir/${sra}_unmapped_read_one -2 $outdir/${sra}_unmapped_read_two -0 $outdir/${sra}_unmapped_read_zero
  #else
  #    printf "File not found!\n" && exit 0
  #fi

  # Determine if SRA is PE or SE. then run magicblast on the nonhost reads. This could be done by looking at the meta data, too FASTA files

  printf "start read count %s\n" "$sra"
  date
  read1_count=$(grep -c "^@" $outdir/${sra}_unmapped.1.fq)
  read2_count=$(grep -c "^@" $outdir/${sra}_unmapped.2.fq)
#  read0_count=$(grep -c "^>" $outdir/${sra}_unmapped_read_zero)

  # Run magicblast using the CARD gene_homology  and CARD variant

  if [[ ( $read2_count > 0 ) && ( $read1_count > 0 ) ]]; then
    printf "Paired End - %d %d\n" "$read1_count" "$read2_count"
    date
    fastx_clipper -i $outdir/${sra}_unmapped.1.fq -o $outdir/${sra}_unmapped_trimmed.1.fq # is this line still needed?
    fastx_clipper -i $outdir/${sra}_unmapped.2.fq -o $outdir/${sra}_unmapped_trimmed.2.fq # is this line still needed?
    printf "Paired End - read1_count %d bowtie align CARDgene sam\n" "$read1_count"
    date
#    magicblast -num_threads $cores  -infmt fasta -query $outdir/${sra}_unmapped_read_one_trimmed -query_mate $outdir/${sra}_unmapped_read_two_trimmed -score 50 -penalty -3 -out $outdir/$sra.CARD_gene.sam -db $cardgene/cardgenedb
#    magicblast -num_threads $cores  -infmt fasta -query $outdir/${sra}_unmapped_read_one_trimmed -query_mate $outdir/${sra}_unmapped_read_two_trimmed -score 50 -penalty -3 -out $outdir/$sra.CARD_snp.sam -db $cardsnp/cardsnpdb
    bowtie2 -x $cardgene/cardgenebt -1 $outdir/${sra}_unmapped_trimmed.1.fq -2 $outdir/${sra}_unmapped_trimmed.2.fq -S $outdir/$sra.CARD_gene_bt.sam
    printf "Paired End - read2_count %d bowtie align CARDsnp sam\n" "$read2_count"
    date
    bowtie2 -x $cardsnp/cardsnpbt -1 $outdir/${sra}_unmapped_trimmed.1.fq -2 $outdir/${sra}_unmapped_trimmed.2.fq -S $outdir/$sra.CARD_snp_bt.sam
#    printf "Paired End - %d magicblast create tabular\n" "$read1_count"
#    date
#    magicblast -num_threads $cores  -infmt fasta -query $outdir/${sra}_unmapped_read_one_trimmed -query_mate $outdir/${sra}_unmapped_read_two_trimmed -score 50 -penalty -3 -outfmt tabular -out $outdir/$sra.CARD_gene.out.tab -db $cardsnp/cardsnpdb
#    magicblast -num_threads $cores  -infmt fasta -query $outdir/${sra}_unmapped_read_one_trimmed -query_mate $outdir/${sra}_unmapped_read_two_trimmed -score 50 -penalty -3 -outfmt tabular -out $outdir/$sra.CARD_snp.out.tab -db $cardsnp/cardsnpdb

  ###JD:in my unpaired run, reads were coming out in read_zero, changed if statement to check in this file rather than read_one, and also clipper and magicblast to read the _zero file
  elif [[ ( $read0_count > 0 ) ]]; then
    printf "Single End - %d %d %d\n" "$read1_count" "$read2_count" "$read0_count"
#    date
#    fastx_clipper -i $outdir/${sra}_unmapped_read_zero -o $outdir/${sra}_unmapped_read_zero_trimmed
    # Run magicblast using CARD Gene Homology nucleotide_fasta_protein_homolog_model.fasta
#    printf "Single End - magicblast CARD gene%d\n" "$read1_count"
#    date
#    magicblast -num_threads $cores -infmt fasta -query $outdir/${sra}_unmapped_read_zero_trimmed -score 50 -penalty -3 -out $outdir/$sra.CARD_gene.sam -db $cardgene/cardgenedb
#    magicblast -num_threads $cores -infmt fasta -query $outdir/${sra}_unmapped_read_zero_trimmed -score 50 -penalty -3 -outfmt tabular -out $outdir/$sra.CARD_gene.out.tab -db $cardgene/cardgenedb
    # Run magicblast using CARD variant  - nucleotide_fasta_protein_variant_model.fasta
    ###JD:corrected cardsnpdb location to $cardsnp/
#    printf "Single End - magicblast CARD SNP%d\n" "$read1_count"
#    date
#    magicblast -num_threads $cores -infmt fasta -query $outdir/${sra}_unmapped_read_zero_trimmed -score 50 -penalty -3 -out $outdir/$sra.CARD_snp.sam -db $cardsnp/cardsnpdb
#    magicblast -num_threads $cores -infmt fasta -query $outdir/${sra}_unmapped_read_zero_trimmed -score 50 -penalty -3 -outfmt tabular -out $outdir/$sra.CARD_snp.out.tab -db $cardsnp/cardsnpdb 
  else
    printf "ERROR no reads to align to CARD - 1 - %d 2 - %d 0 - %d\n" "$read1_count" "$read2_count" "$read0_count" 
  fi

  # Convert the gene_homology aligned SAM to BAM and sort. Delete SAM to save space
  # samtools view -bS $outdir/$sra.CARD_gene.sam | samtools sort -o $outdir/$sra.CARD_gene.bam
  # Convert the gene_variant aligned SAM to BAM and sort. Delete SAM to save space
  # samtools view -bS $outdir/$sra.CARD_snp.sam | samtools sort -o $outdir/$sra.CARD_snp.bam # should we output the unaligned reads to a different bam?
#  rm -f $outdir/$sra.CARD*.sam

  # Optional - Process the above bamfiles into a output-usable form
  printf "Done! %s" "$sra"
  date
  printf "\n"
done > log
