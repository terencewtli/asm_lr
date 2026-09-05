#!/bin/bash
#$ -N W08a_extract_meqtl_cpgs
#$ -cwd
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W08a_extract_meqtl_cpgs.$JOB_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/NA21093/W08a_extract_meqtl_cpgs.$JOB_ID
#
# Extract unique significant meQTL CpG probe IDs from Blueprint mono meQTL file.
# Filters: corrected pval (col 7) < 0.05.
# Output: data/meqtl_mono_sig_cpgs.txt (one probe ID per line)
#
# Blueprint meQTL format (space-delimited, no header):
#   col1: SNP_ID (chr:pos_REF_ALT, hg19)
#   col2: rsID
#   col3: CpG probe ID (450K)
#   col4: nominal pval
#   col5: beta
#   col6: test statistic
#   col7: corrected pval (permutation-based)
#   col8: MAF
#   col9: SE

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr
MEQTL=/u/project/cluo/terencew/claude/project_ideas/me_ca_qtl_seq/data/raw/meqtl_mono.fdr05.tsv.gz
OUTDIR=$PROJDIR/data
mkdir -p $OUTDIR $PROJDIR/logs

echo "Start: $(date)"
echo "Input: $MEQTL"
echo "Extracting CpG probes with corrected pval (col7) < 0.05..."

# stream through 17GB file: filter col7 < 0.05, extract col3, deduplicate
zcat "$MEQTL" 2>/dev/null | \
  awk '$7 < 0.05 {print $3}' | \
  sort -u -T "$OUTDIR" \
  > "$OUTDIR/meqtl_mono_sig_cpgs.txt"

N=$(wc -l < "$OUTDIR/meqtl_mono_sig_cpgs.txt")
echo "Significant meQTL CpGs: $N"

# also save the full significant pairs (SNP + CpG + beta + pval) for effect direction later
echo "Extracting full significant pairs..."
zcat "$MEQTL" 2>/dev/null | \
  awk '$7 < 0.05 {print $1, $3, $4, $5, $7}' | \
  sort -k2,2 -k5,5g -T "$OUTDIR" \
  > "$OUTDIR/meqtl_mono_sig_pairs.txt"

NPAIRS=$(wc -l < "$OUTDIR/meqtl_mono_sig_pairs.txt")
echo "Significant pairs: $NPAIRS"

echo "Done: $(date)"
echo "Outputs:"
echo "  $OUTDIR/meqtl_mono_sig_cpgs.txt   ($N unique CpGs)"
echo "  $OUTDIR/meqtl_mono_sig_pairs.txt  ($NPAIRS pairs)"
