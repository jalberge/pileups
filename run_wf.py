import wolf
import dalmatian
import pandas as pd
import numpy as np

from wolF import *

# extract pairs from a terra workspace
WORKSPACE = "broad-getzlab-mm-germline-t/MM_WGS_GenomeSphere"
wm = dalmatian.WorkspaceManager(WORKSPACE)

WIC = wolf.fc.WorkspaceInputConnector(WORKSPACE)
P = WIC.pairs
S = WIC.get_pairs_as_joint_samples()

PARTICIPANT="PANGEA_4468"

PS = P.merge(S, left_index=True, right_index=True)

PS = PS.loc[ (PS.participant == "PANGEA_4468") & (PS.type_T == "bioskryb")]

# S = S.loc[S.index.str.startswith('Ultra')]

with wolf.Workflow(workflow=forcecall_mafs,
                   scheduler_processes=4,
                   max_concurrent_flows=10,
                   max_concurrent_flow_tasks = 500,
                   common_task_opts = {
                       "retry" : 5,
                       "cleanup_job_workdir" : True
                       }
                   ) as w:
    #for pair, p in PS.iterrows():
    w.run(RUN_NAME="mpileups_test",  # fill in run name
          mafs=PS["mutation_validator_validated_maf_WGS"].tolist(),
          bams=PS["hg38_analysis_ready_bam_T"].tolist(),
          bais=PS["hg38_analysis_ready_bam_index_T"].tolist(),
          samples=PS.index.tolist(),
          n_max=10,
          n_var=2
          )
