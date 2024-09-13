import wolf
import dalmatian
import pandas as pd

from wolF import *


def terra_na(x):
    y = (pd.isna(x)) | (x == "")
    return y


bucket = "gs://acc-genome-sphere/hg38-map-coverage"

bam_table = pd.read_table("fiveprime_bams.tsv", index_col="sample")

with wolf.Workflow(workflow=coverage_interval_list,
                   scheduler_processes=4,
                   max_concurrent_flows=10,
                   max_concurrent_flow_tasks=500,
                   common_task_opts={
                       "retry": 5,
                       "cleanup_job_workdir": False
                   }
                   ) as w:
        w.run(RUN_NAME="cov_haplotype_db",  # fill in run name
              interval_list="https://raw.githubusercontent.com/naumanjaved/fingerprint_maps/master/map_files/hg38_chr.map",
              bams=bam_table["bam"].tolist(),
              bais=bam_table["bai"].tolist(),
              samples=bam_table.index.tolist(),
              localize_bams_to_disks=True,
              n_max=0,
              n_var=250,
              bucket=bucket,
              _debug=True
              )
