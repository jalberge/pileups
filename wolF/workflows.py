import prefect

from .tasks import *


def forcecall_mafs(mafs,
                   bams,
                   bais,
                   samples,
                   n_var=250,
                   n_max=0,
                   fasta="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
                   fasta_index="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
                   fasta_dict="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
                   localize_bams_to_disks=False,
                   bucket=None
                   ):
    @prefect.task
    def sort(x):
        return sorted(x)

    if localize_bams_to_disks:
        local_bams = wolf.LocalizeToDisk(files={"bam": bams, "bai": bais})
        bams = local_bams["bam"]
        bais = local_bams["bai"]
    ref_disk = wolf.LocalizeToDisk(files={"fasta": fasta, "fasta_index": fasta_index, "fasta_dict": fasta_dict})
    variants_lists = Maf2VcfPositions(inputs={
        "mafs": [mafs], "n_var": n_var, "n_max": n_max
    })
    piles = MpileupBams(inputs={
        "bams": [bams],
        "bais": [bais],
        "samples": [samples],
        "fasta": ref_disk["fasta"],
        "fasta_index": ref_disk["fasta_index"],
        "fasta_dict": ref_disk["fasta_dict"],
        "variants_txt": sort(variants_lists["variants"])
    }, overrides = {"bams": "string"} if not localize_bams_to_disks else {})
    ad_dp_matrices = ConcatVcfsToMatrix(inputs={"isec_vcfs": [piles["isec_vcf"]]})
    if bucket is not None:
        wolf.UploadToBucket(
            files=[ad_dp_matrices["tumor_allele_depth"], ad_dp_matrices["total_depth"]],
            bucket=bucket
        )
