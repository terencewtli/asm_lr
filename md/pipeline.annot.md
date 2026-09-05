### Background

User-supplied understanding of the individual steps +
questions/considerations for each step

### Step 1: W00_align_wg.sh

We decided to use individual HG01258 ONT data as the
basis for this project. My understanding is that
HG01258 has three runs of ONT
data, and this script just downloads them, converts to
fastq, then remaps to hg38 and merges the bams.

questions:
- why do we need to remap? are the HPRC default bams
not aligned to a particular reference?
- Why is it named as runs 1/3/5? are there runs 2/4?
Is it normal for individuals to have multiple runs?
- What is the approximate coverage for each run? ONT
data in general? Error rates? Average read length?
Summarize all this for the merged bam?
- Methylation calls are made by guppy? How are these
preserved when the reads are remapped using minimap2?
- More broadly, what does CpG methylation coverage look like
in long read? what does CpG-level coverage look like
genome-wide?

### Step 2: W01_modkit_extract_wg.sh, W02_modkit_pileup_wg.sh

After mapping reads to hg38, we used modkit extract and pileup
to extract CpG level methylation summaries from long-read data.

questions:
- so the modkit extract tsv is every CpG and their read of origin
and methylation status? and modkit pileup just aggregates/summarizes that?
- what is query_kmer in the modkit output?
- how come we used the haplotagged bam for the pileup but not for the
extract step?

### Step 3: W03_phase_wg

After mapping reads, we haplotag them using WhatsHap and the HG12038
reference VCF. that way we get haplotype-specific CpG-level summaries

questions: I think the order should be swapped - we want to run W02
pileup after W03 right?

### Step 4: W05_dasm_loci and W06_merge_chromosomes

After generating per-haplotype pileups, we compare single CpGs and
identify dASM loci. We implement a two-step approach: a fisher test
to identify candidates, and a test that removes loci with methylation
that overlaps between HP1/HP2 (fails to distinguish them)

questions:
- so the basis of the entire results so far is on the nature of the dASM loci. but
just selecting CpGs based on a delta seems a bit fragile. a statistical
test (fisher exact) with multiple testing correction seems to be a more
rigorous way of selecting loci
- are all of the notebooks labeled chr*_W05_dasm_loci.ipynb just the same
thing parallelized across notebooks? does it make sense to wrap this
into a script to run as an array job?

### Step 5:
