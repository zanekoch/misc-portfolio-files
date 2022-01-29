#Ref sequence obtained from https://www.ncbi.nlm.nih.gov/gene/3492 downloaded fasta file 
#contigs obtained from anton's jan 10 run of LJA on hg002
#workflow for aligning reads and calling variants:
# NOTE: manually change reference name to not have spaces otherwise stuff doesn't work
# NOTE: must activate LJA_env conda environment

# hg002_jan10_assembly.fasta
contigs_fn=$1
#  igh_locus_GRCh38.p13_NC_000014.9\:c106879844-105586437.mmi
ref_fn=$2

echo "aligning ${contigs_fn} to ${ref_fn}"
# align
alignment_fn="hg002_jan10_on_igh.sam"

../../../minimap2-2.24_x64-linux/minimap2 -t 20 -a ${2} ${1} > ${alignment_fn}

echo "sorting"
# sort
sorted_alignment_fn="sorted."${alignment_fn}

samtools sort --threads 20 -o $sorted_alignment_fn -O SAM ${alignment_fn}

echo "calling variants"
# call variants and filter
calls_fn="calls.vcf"

bcftools mpileup -f $ref_fn $sorted_alignment_fn | bcftools call -mv -Ov -o $calls_fn

vcfutils.pl varFilter $calls_fn > "filtered."${calls_fn}

# convert to bam
sorted_alignment_bam_fn="${sorted_alignment_fn%.*}.bam"
samtools view --threads 12 -S -b $sorted_alignment_fn > $sorted_alignment_bam_fn

echo "calculating coverage"
# calculate coverage
samtools depth $sorted_alignment_bam_fn > ${sorted_alignment_bam_fn}".coverage"

echo "filtering out reads without large matches"
# filter out reads that do not have atleast 4 digit long matches in cigar string
filtered_sorted_alignment_bam_fn="1k_filtered."${sorted_alignment_bam_fn}
samtools view -h -b $sorted_alignment_bam_fn | awk '($0 ~ /^@/) || ($6 ~ /[0-9]{4,}M/)' > $filtered_sorted_alignment_bam_fn

#re-sort and index these reads
samtools sort --threads 20 -o "sorted."${filtered_sorted_alignment_bam_fn} -O BAM $filtered_sorted_alignment_bam_fn
samtools index "sorted."${filtered_sorted_alignment_bam_fn}
