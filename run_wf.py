import wolf
import dalmatian
import pandas as pd

from wolF import *


def terra_na(x):
    y = (pd.isna(x)) | (x == "")
    return y


bucket = "gs://acc-genome-sphere/hg38-fc-counts"

# extract pairs from a terra workspace
WORKSPACE = "broad-getzlab-mm-germline-t/MM_WGS_GenomeSphere"
wm = dalmatian.WorkspaceManager(WORKSPACE)

WIC = wolf.fc.WorkspaceInputConnector(WORKSPACE)
P = WIC.pairs
S = WIC.get_pairs_as_joint_samples()

PARTICIPANTS = ["PANGEA_4468", "PANGEA_3522", "PANGEA_10704", "PANGEA_3542"]

PS = P.merge(S, left_index=True, right_index=True)

PS = PS.loc[(PS.participant.isin(PARTICIPANTS)) & (PS.type_T == "bioskryb") & ~terra_na(
    PS.mutation_validator_validated_maf_WGS) & ~terra_na(PS.hg38_analysis_ready_bam_T)]

# S = S.loc[S.index.str.startswith('Ultra')]

with wolf.Workflow(workflow=forcecall_mafs,
                   scheduler_processes=4,
                   max_concurrent_flows=10,
                   max_concurrent_flow_tasks=500,
                   common_task_opts={
                       "retry": 5,
                       "cleanup_job_workdir": True
                   }
                   ) as w:
    for PARTICIPANT in PARTICIPANTS:
        rows = PS[PS.participant == PARTICIPANT]
        w.run(RUN_NAME="pileups_maf_" + PARTICIPANT,  # fill in run name
              mafs=rows["mutation_validator_validated_maf_WGS"].tolist(),
              bams=rows["hg38_analysis_ready_bam_T"].tolist(),
              bais=rows["hg38_analysis_ready_bam_index_T"].tolist(),
              samples=rows.index.tolist(),
              localize_bams_to_disks=True,
              n_max=0,
              n_var=250,
              bucket=bucket
              )
