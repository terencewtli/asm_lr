#!/bin/bash
#$ -N W02b_compress_pileups
#$ -l h_data=4G,h_rt=2:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W02b_compress_pileups.$JOB_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr/logs/W02b_compress_pileups.$JOB_ID
#$ -cwd

# bgzip pileup bed files so W05_dasm_loci.ipynb can read them as .bed.gz
# W02 produced .bed files; this step compresses them

PILEUPDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr/pileup
BGZIP=/u/local/apps/htslib/1.12/gcc-4.8.5/bin/bgzip

echo "$(date): compressing pileup beds..."

for bed in $PILEUPDIR/*_hp_pileup_1.bed $PILEUPDIR/*_hp_pileup_2.bed; do
    if [ -f "$bed" ] && [ ! -f "${bed}.gz" ]; then
        echo "  bgzip $bed"
        $BGZIP -@ 4 -k "$bed"
    fi
done

echo "$(date): done"
ls -lh $PILEUPDIR/*.bed.gz 2>/dev/null | wc -l
echo "bed.gz files created"
