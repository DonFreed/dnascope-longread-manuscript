#!/usr/bin/env snakemake

"""
Data processing pipeline for the DNAscope LongRead paper
"""

# Copyright (c) 2022 Sentieon Inc. All rights reserved

CONFIG_PATH = f"{workflow.snakefile}.yaml"
configfile: CONFIG_PATH

import re

FEMALE_SEX = "female"
MALE_SEX = "male"

###### Output files ######
## Data files ##
DIPLOID_BED = "diploid_bed/{ref}/{sex}/diploid.bed"

## Preprocessing ##
PBMM2 = "pbmm2/{ref}/{sample}/{aln_params}/{read_idx}/aligned.bam"
MERGE = "merged_aln/{ref}/{sample}/{aln_params}/merged.cram"
SUBSET = "subset/{ref}/{sample}/{aln_params}/{subset}/subset.cram"

## Metrics ##
COVERAGE = "coverage/{ref}/{sample}/{aln_params}/{subset}/coverage"

## Calling and Evaluation ##
DNASCOPE = "dnascope/{ref}/{sample}/{aln_params}/{subset}/calls.vcf.gz"
HAPPY = "happy/{ref}/{sample}/{aln_params}/{subset}/{truthset}/eval"
MERGED_HAPPY = "happy_merged/happy_extended_metrics.tsv"

## Figures and tables ##
FIGURES = "analysis/figures"


###### Filename patterns #####
happy_extended_csv_pat = re.compile(
    "happy/(?P<ref>.*)/(?P<sample>.*)/(?P<aln_params>.*)/(?P<subset>.*)/(?P<truthset>.*)/eval.extended.csv"
)

###### Input/Parameter functions ######
def get_all_happy_evals(suffix=".summary.csv"):
  expected = []

  # Static wildcards
  ref = "hs38"
  aln_params = "--preset HiFi -c 0 -y 70"

  # All samples - full depth
  for sample in config["input"]["samples"].keys():
    for subset in ("full",):
      truthset = "v4.2.1"
      expected.append(HAPPY.format(**locals()) + suffix)
    
  # pFDA HG003 - serial downsamples
  for sample in ("pFDA-truth_V2-PB-HG003",):
    # Subsets from 5x to 30x with a 5x step
    for subset in config["pipeline"]["hg003_subsets"]:
      truthset = "v4.2.1"
      expected.append(HAPPY.format(**locals()) + suffix)

  # Evaluate the Google Health pFDA Truth V2 submissions against the latest truthset
  subset = "full"
  truthset = "v4.2.1"
  for truth_sample in ("HG002", "HG003", "HG004"):
    sample = f"GH_pFDA_truth_V2--{truth_sample}"
    expected.append(HAPPY.format(**locals()) + suffix)

  # pFDA HG002 - CMRG benchmark
  sample = "pFDA-truth_V2-PB-HG002"
  subset = "full"
  truthset = "CMRGv1.00"
  expected.append(HAPPY.format(**locals()) + suffix)

  return expected


def get_all_coverage_metrics():
  expected = []
  
  ref = "hs38"
  aln_params = "--preset HiFi -c 0 -y 70"

  # All samples - full depth
  subset = "full"
  for sample in config["input"]["samples"].keys():
    expected.append(COVERAGE.format(**locals()) + ".sample_summary")

  # pFDA HG003 - serial downsamples
  sample = "pFDA-truth_V2-PB-HG003"
  for subset in config["pipeline"]["hg003_subsets"]:
    expected.append(COVERAGE.format(**locals()) + ".sample_summary")

  return expected


def get_pbmm2_reads(wildcards):
  sample_d = config["input"]["samples"][wildcards.sample]
  idx = int(wildcards.read_idx)
  if "fq" in sample_d:
    return sample_d["fq"][idx]
  elif "bam" in sample_d:
    return sample_d["bam"][idx]
  else:
    raise ValueError(f"Cannot find reads from sample dict: {sample_d}")


def get_pbmm2_rg(wildcards):
  sample_d = config["input"]["samples"][wildcards.sample]
  idx = int(wildcards.read_idx)
  if "rg" in sample_d:  # Sample is a 'fq' in the config
    return "--rg '" + sample_d["rg"][idx] + "'"
  else:
    return ""


def get_merge_alns(suffix=""):
  def _get_merge_alns(wildcards):
    sample_d = config["input"]["samples"][wildcards.sample]
    n_alns = -1
    if "rg" in sample_d:
      n_alns = len(sample_d["rg"])
    elif "bam" in sample_d:
      n_alns = len(sample_d["bam"])
    else:
      raise ValueError("Unknown number of alignment files")

    pbmm2_alns = []
    for read_idx in range(n_alns):
      pbmm2_alns.append(PBMM2.format(read_idx=read_idx, **wildcards.__dict__) + suffix)
    return pbmm2_alns
  return _get_merge_alns


def get_dnascope_diploid_bed(wildcards):
  sample = config["sample_map"][wildcards.sample]
  gender = config["gender_map"][sample]
  return DIPLOID_BED.format(ref=wildcards.ref, sex=gender)


def get_calls_vcf(suffix=""):
  def _get_calls_vcf(wldc):
    if "GH_pFDA_truth_V2" in wldc.sample:
      for truth_sample in ("HG002", "HG003", "HG004"):
        if truth_sample in wldc.sample:
          return config["input"]["pFDA_Truth_V2_submissions"]["Google_Health"][truth_sample] + suffix
      else:
        raise ValueError(f"Cannot determine truth sample from: {wldc.sample}")
    return DNASCOPE.format(**wldc.__dict__) + suffix
  return _get_calls_vcf


def get_truth_vcf(suffix=""):
  def _get_truth_vcf(wldc):
    truth_sample = config["sample_map"][wldc.sample]
    truth_ref = config["ref_map"][wldc.ref]
    return config["input"]["truth"][truth_ref][wldc.truthset][truth_sample]["vcf"] + suffix
  return _get_truth_vcf


def get_truth_bed(wldc):
  truth_sample = config["sample_map"][wldc.sample]
  truth_ref = config["ref_map"][wldc.ref]
  return config["input"]["truth"][truth_ref][wldc.truthset][truth_sample]["bed"]


def get_collect_sample_names(wldc, input):
  """ Build coherent sample names from the wildcards of the input files """
  sample_names = []
  for extended_csv in input.happy_evals:
    m = happy_extended_csv_pat.match(extended_csv)
    if not m:
      raise ValueError(f"Extended csv output file, '{extended_csv}' did not match the expected pattern")
    wildcards = m.groupdict()
    # Build the sample name from the variable wildcards
    sample_name = f"{wildcards['sample']}--{wildcards['subset']}--{wildcards['truthset']}"
    sample_names.append(sample_name)

  # Separate the sample names for the shell
  sample_names = "'" + "' '".join(sample_names) + "'"
  return sample_names


###### Rules ######
rule all:
  input:
    FIGURES + ".html",

# Generate a BED file of diploid regions in the genome
rule diploid_bed:
  input:
    ref_fai = lambda wldc: config["input"]["ref"][wldc.ref] + ".fai",
  output:
    bed = DIPLOID_BED,
  log:
    DIPLOID_BED + ".log",
  threads:
    8
  params:
    extra_chr = lambda wldc: "chrX" if wldc.sex == FEMALE_SEX else ""
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.bed}")
    mkdir -p "$outdir"
    exec &>"{log}"

    # Get the list of chromosomes. Only hg19/hg38 style 'chr' is supported
    tmp_fa_file="${{outdir}}"/tmp.txt
    rm "$tmp_fa_file" || $(exit 0)
    for i in $(seq 1 22); do
      echo "^chr$i	" >> "$tmp_fa_file"
    done
    if [[ "{params.extra_chr}" == "chrX" ]]; then
      echo "^chrX	" >> "$tmp_fa_file"
    fi

    grep -f "$tmp_fa_file" "{input.ref_fai}" | \
      awk -v 'OFS=\t' '{{print $1,0,$2}}' > "{output.bed}"
    rm "$tmp_fa_file"
    """

# Align reads with pbmm2
rule pbmm2:
  input:
    reads = get_pbmm2_reads,
    ref = lambda wldc: config["input"]["ref"][wldc.ref],
    pbmm2 = config["tools"]["pbmm2"],
  output:
    aln = temp(PBMM2),
    idx = temp(PBMM2 + ".bai"),
  log:
    stdout = PBMM2 + ".stdout",
    stderr = PBMM2 + ".stderr",
  threads:
    32
  params:
    rg = get_pbmm2_rg,
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.aln}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}"

    {input.pbmm2} align --sort -j {threads} {params.rg} \
      {wildcards.aln_params} "{input.ref}" "{input.reads}" "{output.aln}"
    """

# Merge alignments
rule merge_alns:
  input:
    aln = get_merge_alns(),
    idx = get_merge_alns(".bai"),
    ref = lambda wldc: config["input"]["ref"][wldc.ref],
    sentieon = config["tools"]["sentieon"],
  output:
    aln = MERGE,
    idx = MERGE + ".crai",
  log:
    stdout = MERGE + ".stdout",
    stderr = MERGE + ".stderr",
  threads:
    32
  params:
    input_alns = lambda wldc, input: " -i '" + "' -i '".join(input.aln) + "'",
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.aln}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    {input.sentieon} driver {params.input_alns} -t {threads} -r {input.ref} \
      --algo ReadWriter --cram_write_options version=3.0,compressor=rans \
      "{output.aln}"
    """

rule subset_alns:
  input:
    aln = rules.merge_alns.output.aln,
    idx = rules.merge_alns.output.idx,
    ref = lambda wldc: config["input"]["ref"][wldc.ref],
    samtools = config["tools"]["samtools"],
  output:
    aln = SUBSET,
    idx = SUBSET + ".crai",
  log:
    stdout = SUBSET + ".stdout",
    stderr = SUBSET + ".stderr",
  threads:
    16
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.aln}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    if [[ "{wildcards.subset}" = "full" ]]; then
      # No subset - hard link input to output
      ln "{input.aln}" "{output.aln}"
      ln "{input.idx}" "{output.idx}"
    else
      # Take a random subset of the reads
      random_seed=$RANDOM
      {input.samtools} view --subsample-seed $random_seed \
        --subsample 0.{wildcards.subset} -O CRAM -T {input.ref} \
        -o "{output.aln}" -@ {threads} --output-fmt-option version=3.0 "{input.aln}"
      {input.samtools} index -c "{output.aln}"
    fi
    """

rule coverage:
  input:
    aln = rules.subset_alns.output.aln,
    idx = rules.subset_alns.output.idx,
    ref = lambda wldc: config["input"]["ref"][wldc.ref],
    sentieon = config["tools"]["sentieon"],
  output:
    sample_summary = COVERAGE + ".sample_summary",
  log:
    COVERAGE + ".log",
  threads:
    16
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.sample_summary}")
    mkdir -p "$outdir"
    exec &>"{log}"

    {input.sentieon} driver -t {threads} -r {input.ref} -i "{input.aln}" \
      --algo CoverageMetrics --omit_base_output --omit_locus_stat \
      --omit_interval_stat "$outdir"/coverage
    """


rule dnascope_hifi:
  input:
    aln = rules.subset_alns.output.aln,
    idx = rules.subset_alns.output.idx,
    ref = lambda wldc: config["input"]["ref"][wldc.ref],
    sentieon = config["tools"]["sentieon"],
    bedtools = config["tools"]["bedtools"],
    bcftools = config["tools"]["bcftools"],
    dnascope_hifi = config["tools"]["dnascope_hifi"],
    dnascope_hifi_model = config["tools"]["dnascope_hifi_model"],
    diploid_bed = get_dnascope_diploid_bed,
    mhc_bed = lambda wldc: config["input"]["mhc"][wldc.ref],
    dbsnp_vcf = lambda wldc: config["input"]["dbsnp"][wldc.ref],
  output:
    vcf = DNASCOPE,
    tbi = DNASCOPE + ".tbi",
  benchmark:
    DNASCOPE + ".benchmark.txt",
  log:
    stdout = DNASCOPE + ".stdout",
    stderr = DNASCOPE + ".stderr",
  threads:
    64
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.vcf}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    # Move the input BAM onto the local SSD for an accurate runtime benchmark
    tmpdir=$(mktemp -d -p /data/work/dfreed/tmp)
    cp "{input.aln}" "$tmpdir"/sample.cram
    cp "{input.aln}".crai "$tmpdir"/sample.cram.crai

    # Set the PATH correctly for the required tools
    sentieon_dir=$(dirname "{input.sentieon}")
    bcftools_dir=$(dirname "{input.bcftools}")
    bedtools_dir=$(dirname "{input.bedtools}")
    export PATH="$sentieon_dir":"$bcftools_dir":"$bedtools_dir":$PATH
    echo "Start time: $(date +'%s')"
    bash "{input.dnascope_hifi}" -d "{input.dbsnp_vcf}" -B "{input.mhc_bed}" \
      -m "{input.dnascope_hifi_model}" -b "{input.diploid_bed}" -t {threads} \
      -r {input.ref} -i "$tmpdir"/sample.cram "$tmpdir"/calls.vcf.gz
    echo "End time: $(date +'%s')"

    # Move the output VCF and remove the tmp directory
    mv "$tmpdir"/calls.vcf.gz "{output.vcf}"
    mv "$tmpdir"/calls.vcf.gz.tbi "{output.tbi}"
    rm -r "$tmpdir"
    """

rule happy:
  input:
    calls_vcf = get_calls_vcf(),
    calls_tbi = get_calls_vcf(".tbi"),
    truth_vcf = get_truth_vcf(),
    truth_tbi = get_truth_vcf(".tbi"),
    truth_bed = get_truth_bed,
    ref = lambda wldc: config["input"]["ref"][wldc.ref],
    sdf = lambda wldc: config["input"]["sdf"][wldc.ref],
    strat_tsv = lambda wldc: config["input"]["strat_tsv"][config["ref_map"][wldc.ref]],
    python2 = config["tools"]["python2"],
    happy = config["tools"]["happy"],
    rtg_dir = config["tools"]["rtg_dir"],
  output:
    summary = HAPPY + ".summary.csv",
    extended = HAPPY + ".extended.csv",
    vcf = HAPPY + ".vcf.gz",
    tbi = HAPPY + ".vcf.gz.tbi",
  log:
    stdout = HAPPY + ".stdout",
    stderr = HAPPY + ".stderr",
  threads:
    8
  resources:
    mem_mb=36000,
  params:
    gender = lambda wldc: config["gender_map"][config["sample_map"][wldc.sample]],
    strat_str = lambda wldc, input: "--stratification " + input.strat_tsv,
    happy_xargs = "--no-decompose --no-leftshift",
  shell:
    """
    set -exvo pipefail
    outdir=$(dirname "{output.summary}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}"

    # Link the input file into the output directory and cd to outdir. Hap.py 
    # cannot handle complex path characters (', =, ' ', etc.)
    in_vcf=$(realpath "{input.calls_vcf}")
    rm "$outdir"/calls.vcf.gz "$outdir"/calls.vcf.gz.tbi || $(exit 0)
    ln -s "$in_vcf" "$outdir"/calls.vcf.gz
    ln -s "$in_vcf".tbi "$outdir"/calls.vcf.gz.tbi

    outpre=$(basename "{output.summary}")
    outpre="${{outpre%%.summary.csv}}"

    cd "$outdir"

    export PATH={input.rtg_dir}:$PATH
    {input.python2} {input.happy} {input.truth_vcf} \
      calls.vcf.gz --verbose -o "$outpre" -r {input.ref} --gender \
      {params.gender} -f {input.truth_bed} --threads \
      {threads} --engine=vcfeval --engine-vcfeval-template {input.sdf} \
      {params.strat_str} {params.happy_xargs}
    """

rule merge_happy:
  input:
    happy_evals = get_all_happy_evals(".extended.csv"),
    python3 = config["tools"]["python3"],
    extract_extended_metrics = config["tools"]["extract_extended_metrics"],
  output:
    MERGED_HAPPY,
  log:
    MERGED_HAPPY + ".log",
  threads:
    4
  params:
    infiles = lambda wldc, input: " '" + "' '".join(input.happy_evals) + "'",
    sample_names = get_collect_sample_names,
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output}")
    mkdir -p "$outdir"
    exec &>"{log}"

    {input.python3} {input.extract_extended_metrics} \
      --sample_names {params.sample_names} \
      --infiles {params.infiles} \
      --filters "PASS" \
      --subtype '*' \
      --columns "Type" "Subset" "METRIC.Recall" "METRIC.Precision" \
        "METRIC.F1_Score" "TRUTH.TOTAL" "TRUTH.TP" "TRUTH.FN" "QUERY.FP" \
        "TOTAL.ERRORS" \
      > "{output}"
    """

rule plot_figures:
  input:
    happy_evals = rules.merge_happy.output[0],
    coverage_metrics = get_all_coverage_metrics(),
    figure_notebook = config["tools"]["figure_notebook"],
    papermill = config["tools"]["papermill"],
    jupyter = config["tools"]["jupyter"],
  output:
    ipynb = FIGURES + ".ipynb",
    html = FIGURES + ".html",
  params:
    config_path = CONFIG_PATH,
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output}")
    mkdir -p "$outdir"
    
    {input.papermill} {input.figure_notebook} "{output.ipynb}" --prepare-only \
      -p 'work_dir' "$(pwd)" \
      -p 'dnascope_output_format' 'dnascope/{{ref}}/{{sample}}/{{aln_params}}/{{subset}}/calls.vcf.gz' \
      -p 'happy_merged_output' 'happy_merged/happy_extended_metrics.tsv' \
      -p 'config_file' '{params.config_path}'

    {input.jupyter} nbconvert --execute --to html "{output.ipynb}" \
      --ExecutePreprocessor.timeout=-1
    """
