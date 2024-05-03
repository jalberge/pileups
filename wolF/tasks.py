import wolf

samtools_docker = "gcr.io/broad-getzlab-workflows/samtools@sha256" \
                  ":8074df347e20ca7f39646914eb9fcccde30a9ca85d4dd10b60fe99d5b78a223d "


class Maf2VcfPositions(wolf.Task):
    name = "Maf2VcfPositions"
    inputs = {"mafs": None, "n_var": 250, "n_max": 0}
    script = """

set -euxo pipefail

# use n_max>0 for test purposes

# SBS only for now
# mafs is a list of maf files!
# first extract positions; sort; dedup; split

xargs cut -f5,6,11,13 < $mafs | \
    awk -v OFS="\t" '$3 ~ /^[ACGT]$/ && $4 ~ /^[ACGT]$/' | \
    sort -k1,2 -V | \
    uniq > master

if [[ $n_max -gt 0 ]]; then shuf -n $n_max master | sort -k1,2 -V > master_ && mv master_ master; fi

split -l $n_var -a 5 -d --additional-suffix=.txt master variants_
    
    """
    outputs = {
        "variants": "variants_*.txt",
    }
    docker = samtools_docker


class MpileupBams(wolf.Task):
    name = "MpileupBams"
    inputs = {
        "bams": None,
        "bais": None,
        "samples": None,
        "fasta": None,
        "fasta_index": None,
        "fasta_dict": None,
        "variants_txt": None
    }
    # overrides = {"bams": "string"}
    script = """
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

set -euxo pipefail

# ln -s indexes from common > inputs > workspace
for x in $(cat $bais); do ln -s $x . ; done

shard=$(basename $variants_txt .txt)

# prepare dummy vcf header
# also extract the positions for mpileup / view streams
echo -e "##fileformat=VCFv4.3" > variants.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> variants.vcf
awk -v OFS="\t" '{print $1,$2,".",$3,$4,".",".","."}' $variants_txt >> variants.vcf
bgzip variants.vcf
tabix -p vcf variants.vcf.gz

# extract positions to stream bam
cut -f1,2 $variants_txt > positions.txt

# convert file to list of bams
bam_list=$(cat $bams | tr "\n" " ")

# create sample map
paste -d" " $bams $samples > sample_map

# mpileup depth and t_alt_count
# norm to split multi allele
bcftools mpileup -a FORMAT/AD,FORMAT/DP -A -d 100 -R positions.txt --ignore-RG -I -f $fasta $bam_list | \
    bcftools reheader -s sample_map | \
    bcftools norm -m - --write-index -o bcfpiles.vcf.gz
# isec to intersect with master list of variants (keep only ALT allele and exclude *)
bcftools isec -c none -p . -n=2 -w1  bcfpiles.vcf.gz variants.vcf.gz

mv 0000.vcf ${shard}_isec.vcf

# clean indexes asap
# TODO find a way to link the indexes
rm *.bai
"""
    outputs = {
        "isec_vcf": "*_isec.vcf"
    }
    docker = samtools_docker
    # resources = {"mem": "8G"}  # mpileup 1.20 is not really parallelized. just compression vcf.


# define additional tasks in the same way that task1 is defined above.

class ConcatVcfsToMatrix(wolf.Task):
    name = "ConcatVcfsToMatrix"
    inputs = {
        "isec_vcfs": None
    }
    script = """
    set -exuo pipefail

bcftools concat -f $isec_vcfs -o concat_0000.vcf

    # query to transform VCF to MTX
bcftools query -f '%CHROM\_%POS\_%REF\_%ALT[\t%AD{1}]' concat_0000.vcf > AD.txt
bcftools query -f '%CHROM\_%POS\_%REF\_%ALT[\t%DP]' concat_0000.vcf > DP.txt
bcftools query -l concat_0000.vcf > samples.txt

# rm last tab created by pivot
cat samples.txt | tr "\n" "\t" | sed 's/\t$//' > header

awk '1' header DP.txt  > final_dp.txt
awk '1' header AD.txt  > final_ad.txt
"""
    outputs = {
        "samples": "samples.txt",
        "tumor_allele_depth": "final_ad.txt",
        "total_depth": "final_dp.txt"
    }
    docker = samtools_docker
