#!/bin/bash
#$ -N W05_merge_wg
#$ -l h_data=8G,h_rt=1:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W05_merge_wg.$JOB_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W05_merge_wg.$JOB_ID
#$ -hold_jid W04_dasm_loci_wg
#$ -cwd

# Merge per-chromosome dASM tables → whole-genome parquet,
# then execute the W05_dasm_loci summary notebook (single output).
# Runs after W04 job array completes.

source /u/home/t/terencew/project-cluo/miniconda3/bin/activate allcools

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
SAMPLE=NA21093

echo "$(date): W05 — merging per-chromosome tables"
python3 - <<'PYEOF'
import glob, json, os
import pandas as pd

projdir = '/u/project/cluo/terencew/claude/project_ideas/asm_lr'
paths = sorted(glob.glob(f'{projdir}/tables/phasing_results_chr*.parquet'))
if not paths:
    raise FileNotFoundError('No phasing_results_chr*.parquet files found')

missing = [p for p in paths if not os.path.exists(p)]
if missing:
    print(f'WARNING: missing {len(missing)} files: {missing}')

dfs = [pd.read_parquet(p) for p in paths]
wg  = pd.concat(dfs, ignore_index=True)
wg.to_parquet(f'{projdir}/tables/dasm_loci_wg.parquet', index=False)

n_phased = int(wg['read_level_support'].sum())
summary = {
    'sample':        'NA21093',
    'n_dasm_loci':   int(len(wg)),
    'n_phased':      n_phased,
    'pct_phased':    round(n_phased / len(wg) * 100, 2),
    'n_in_rosenski': int(wg['in_rosenski'].sum()),
}
with open(f'{projdir}/tables/phasing_summary_wg.json', 'w') as f:
    json.dump(summary, f, indent=2)
print(f'Merged {len(wg):,} loci ({n_phased:,} phased) → dasm_loci_wg.parquet')
import json as _j; print(_j.dumps(summary, indent=2))
PYEOF

echo "$(date): running W05 summary notebook"
NB_IN=$PROJDIR/notebooks/wg/W05_dasm_loci.ipynb
NB_OUT=$PROJDIR/notebooks/wg/executed/W05_dasm_loci.ipynb
mkdir -p $PROJDIR/notebooks/wg/executed

papermill \
    --kernel allcools \
    --log-output \
    -p projdir $PROJDIR \
    -p sample  $SAMPLE \
    $NB_IN \
    $NB_OUT

echo "$(date): W05 done — see notebooks/wg/executed/W05_dasm_loci.ipynb"
