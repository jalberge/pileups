import prefect

from .tasks import *


def forcecall_mafs(mafs,
                   bams,
                   bais,
                   samples,
                   participant="P",
                   n_var=250,
                   n_max=0,
                   fasta="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
                   fasta_index="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
                   fasta_dict="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
                   localize_bams_to_disks=False,
                   bucket=None,
                   project="dfci-ghobriallab-gcp",
                   workspace=None,
                   workspace_entity_name=None,
                   workspace_entity_type="participant",
                   _debug=False
                   ):
    @prefect.task
    def sort(x):
        return sorted(x)

    @prefect.task
    def atleast1d(x):
        return [x] if (not isinstance(x, list)) and x is not None else x

    if localize_bams_to_disks:
        local_bams = wolf.LocalizeToDisk(files={"bam": bams, "bai": bais}, project=project)
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
        "variants_txt": sort(atleast1d(variants_lists["variants"])),
        "isec_collapse_mode": "none"  # exact ALT match
    }, overrides={"bams": "string"} if not localize_bams_to_disks else {})

    ad_dp_matrices = ConcatVcfsToMatrix(inputs={"isec_vcfs": [piles["isec_vcf"]]})

    if localize_bams_to_disks and not _debug:
        wolf.localization.DeleteDisk(
            name="DeleteBams",
            inputs={
                "disk": local_bams["bam"],
                "upstream": ad_dp_matrices["total_depth"]
            },
            mapped=True
        )

    if bucket is not None:
        wolf.UploadToBucket(
            files=[ad_dp_matrices["tumor_allele_depth"], ad_dp_matrices["total_depth"], ad_dp_matrices["samples"]],
            bucket=bucket+'/'+str(participant)
        )

    if workspace is not None:
        if workspace_entity_name is None:
            workspace_entity_name = participant
        attr_map = {
                "mut_pileups_tumor_allele_depth" : ad_dp_matrices["tumor_allele_depth"],
                "mut_pileups_total_depth" : ad_dp_matrices["total_depth"],
                "mut_pileups_samples" : ad_dp_matrices["samples"]
                }
        wolf.fc.SyncToWorkspace(
                nameworkspace = workspace,
                entity_type = workspace_entity_type,
                entity_name = workspace_entity_name,
                attr_map = attr_map
        )

def coverage_interval_list(interval_list,
                           bams,
                           bais,
                           samples,
                           n_var=250,
                           n_max=0,
                           fasta="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
                           fasta_index="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
                           fasta_dict="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
                           localize_bams_to_disks=False,
                           bucket=None,
                           _debug=False
                           ):
    @prefect.task
    def sort(x):
        return sorted(x)

    @prefect.task
    def atleast1d(x):
        return [x] if (not isinstance(x, list)) and x is not None else x

    if localize_bams_to_disks:
        local_bams = wolf.LocalizeToDisk(files={"bam": bams, "bai": bais})
        bams = local_bams["bam"]
        bais = local_bams["bai"]

    ref_disk = wolf.LocalizeToDisk(files={"fasta": fasta, "fasta_index": fasta_index, "fasta_dict": fasta_dict})

    variants_lists = IntervalList2VcfPositions(inputs={
        "interval_list": [interval_list], "n_var": n_var, "n_max": n_max
    })

    piles = MpileupBams(inputs={
        "bams": [bams],
        "bais": [bais],
        "samples": [samples],
        "fasta": ref_disk["fasta"],
        "fasta_index": ref_disk["fasta_index"],
        "fasta_dict": ref_disk["fasta_dict"],
        "variants_txt": sort(atleast1d(variants_lists["variants"])),
        "isec_collapse_mode": "all"  # don't need ALT match
    }, overrides={"bams": "string"} if not localize_bams_to_disks else {})

    ad_dp_matrices = ConcatVcfsToMatrix(inputs={"isec_vcfs": [piles["isec_vcf"]]})

    if localize_bams_to_disks and not _debug:
        wolf.localization.DeleteDisk(
            name="DeleteBams",
            inputs={
                "disk": local_bams["bam"],
                "upstream": ad_dp_matrices["total_depth"]
            },
            mapped=True
        )

    if bucket is not None:
        wolf.UploadToBucket(
            files=[ad_dp_matrices["tumor_allele_depth"], ad_dp_matrices["total_depth"], ad_dp_matrices["samples"]],
            bucket=bucket
        )

def genotype_in_the_cloud(bam, bai, name,

                          fasta="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
                          fasta_index="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
                          fasta_dict="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",

                          bucket="gs://acc-genome-sphere/hg38-5prime-gt/",
                          version="hg38",
                          regions="gs://jba-utils/hg38_5prime_20240710_n_1913_positions_sorted.txt",

                          workspace = None,
                          workspace_entity_name=None,
                          workspace_entity_type="sample"):

    ref_disk = wolf.LocalizeToDisk(files={"fasta": fasta, "fasta_index": fasta_index, "fasta_dict": fasta_dict})

    genotypes = BcftoolsMpileupCallCloud(inputs={

        "bam":bam,"bai":bai,

        "fasta":ref_disk["fasta"],
        "fasta_index":ref_disk["fasta_index"],
        "fasta_dict":ref_disk["fasta_dict"],

        "regions": regions,
        "version":version,

        "name":name,
    })

    if bucket is not None:
        wolf.UploadToBucket(
            files=[ genotypes["vcf_gz"], genotypes["vcf_gz_tbi"]],
            bucket=bucket
        )

    if workspace is not None:
        if workspace_entity_name is None:
            workspace_entity_name = name
        attr_map = {
                "fingerprints_fiveprime_vcf_gz" : genotypes["vcf_gz"],
                "fingerprints_fiveprime_vcf_gz_tbi" : genotypes["vcf_gz_tbi"],
                }
        wolf.fc.SyncToWorkspace(
                nameworkspace = workspace,
                entity_type = workspace_entity_type,
                entity_name = workspace_entity_name,
                attr_map = attr_map
        )

def merge_genotypes(vcf_gz_list, vcf_gz_tbi_list, sample_set, bucket="gs://acc-genome-sphere/hg38-5prime-gt-merge/",
                    workspace=None, workspace_entity_name=None, workspace_entity_type="sample_set"):

    merged_genotypes = MergeVcfs(inputs={
        "vcf_gz_list": [vcf_gz_list],
        "vcf_gz_tbi_list": [vcf_gz_tbi_list],
        "sample_set": sample_set
    })

    if bucket is not None:
        wolf.UploadToBucket(
            files=[ genotypes["vcf_gz"], genotypes["vcf_gz_index"]],
            bucket=bucket
        )

    if workspace is not None:
        if workspace_entity_name is None:
            workspace_entity_name = sample_set
        attr_map = {
            "merged_genotypes_vcf_gz": merged_genotypes["vcf_gz"],
            "merged_genotypes_vcf_gz_tbi": merged_genotypes["vcf_gz_tbi"],
        }
        wolf.fc.SyncToWorkspace(
            nameworkspace=workspace,
            entity_type=workspace_entity_type,
            entity_name=workspace_entity_name,
            attr_map=attr_map
        )
