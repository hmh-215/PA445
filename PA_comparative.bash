#!/bin/bash
# P. aeruginosa comparative genomics - PA445 vs global blaNDM-1 strains
# Workflow:
#   Phase 1  : Download full dataset (NDM + background), MLST/Pasty/RGI/Prokka/Panaroo
#   Phase 2  : prokka_rgi_reconcile.py -> identify true NDM-1 strains
#   Phase 3  : Subsample dataset (all NDM-1 + =5 background per country)
#   Phase 4  : RGI heatmap (subsampled) + IQ-TREE phylogeny (subsampled)
#   Phase 5  : gbk_extraction + clinker on NDM-1 strains only
#   Phase 6  : Plasmid analysis on NDM-1 strains only
#   Phase 7  : Subsubsampling dataset analyses
#   Phage 8  : Building time-caliberated phylogenetic tree with BEAST
# Building no. 4
set -euo pipefail

################
# GLOBAL SETUP #
################

SAMPLE_PATH="/storage/student9/projects/Paeruginosa"
REF_PATH="/storage/student9/references"
COMP_PATH="${SAMPLE_PATH}/comparative_genomics_Paeruginosa"
REFERENCES_PA="${COMP_PATH}/references"
TOOL_PATH="/storage/student9/tools"

MLST_PA="${COMP_PATH}/mlst"
PASTY_PA="${COMP_PATH}/pasty"
ANI_PA="${COMP_PATH}/ani"
RGI_PA="${COMP_PATH}/rgi"
PROKKA_PA="${COMP_PATH}/prokka"
PANAROO_PA="${COMP_PATH}/panaroo"
ANI_PA="${COMP_PATH}/ani"
IQTREE_PA="${COMP_PATH}/iqtree"
GENE_PA="${COMP_PATH}/gene"
CLINKER_PA="${COMP_PATH}/clinker"
SUBSAMPLE_PA="${COMP_PATH}/subsampled"
SUBSAMPLE_FASTA="${SUBSAMPLE_PA}/fasta"
RGI_SUB="${COMP_PATH}/rgi_subsampled"
IQTREE_SUB="${COMP_PATH}/iqtree_subsampled"
PLASMID_PA="${COMP_PATH}/plasmid"
PLATON_PA="${PLASMID_PA}/platon"
PLASMIDFINDER_PA="${PLASMID_PA}/plasmidfinder"
MOBSUITE_PA="${PLASMID_PA}/mobsuite"
SUBSUB_PA="${COMP_PATH}/subsubsampled"
SUBSUB_FASTA="${SUBSUB_PA}/fasta"
RGI_SUBSUB="${COMP_PATH}/rgi_subsubsampled"
PANAROO_SUBSUB="${COMP_PATH}/panaroo_subsubsampled"
IQTREE_SUBSUB="${COMP_PATH}/iqtree_subsubsampled"
BEAST_PA="${COMP_PATH}/beast_dating"
SNIPPY_PA="${BEAST_PA}/snippy"
GUBBINS_PA="${BEAST_PA}/gubbins"

mkdir -p "${COMP_PATH}/references"
mkdir -p "${COMP_PATH}/mlst"
mkdir -p "${COMP_PATH}/pasty"
mkdir -p "${COMP_PATH}/ani"
mkdir -p "${COMP_PATH}/rgi"
mkdir -p "${RGI_PA}/rgi_inputs"
mkdir -p "${COMP_PATH}/prokka"
mkdir -p "${COMP_PATH}/panaroo"
mkdir -p "${ANI_PA}/full"
mkdir -p "${ANI_PA}/subsampled"
mkdir -p "${ANI_PA}/subsubsampled"
mkdir -p "${PANAROO_PA}/gff_inputs"
mkdir -p "${COMP_PATH}/gene"
mkdir -p "${COMP_PATH}/clinker"
mkdir -p "${CLINKER_PA}/truepos_gbk"
mkdir -p "${CLINKER_PA}/gbk_inputs"
mkdir -p "${COMP_PATH}/plasmid"
mkdir -p "${PLASMID_PA}/platon"
mkdir -p "${PLASMID_PA}/plasmidfinder"
mkdir -p "${PLASMID_PA}/mobsuite"
mkdir -p "${SUBSAMPLE_PA}"
mkdir -p "${SUBSAMPLE_FASTA}"
mkdir -p "${RGI_SUB}/rgi_inputs"
mkdir -p "${IQTREE_SUB}"
mkdir -p "${REFERENCES_PA}/download_logs"
mkdir -p "${COMP_PATH}/subsubsampled"
mkdir -p "${SUBSUB_PA}/fasta"
mkdir -p "${COMP_PATH}/rgi_subsubsampled"
mkdir -p "${RGI_SUBSUB}/rgi_inputs"
mkdir -p "${PANAROO_SUBSUB}/results"
mkdir -p "${PANAROO_SUBSUB}/gff_inputs"
mkdir -p "${COMP_PATH}/iqtree_subsubsampled"
mkdir -p "${COMP_PATH}/beast_dating"
mkdir -p "${BEAST_PA}/snippy"
mkdir -p "${BEAST_PA}/gubbins"

#Set variables
PA445="${SAMPLE_PATH}/445.hybrid_assembly_filter500bp.renamed.fasta"
REF_ANCHORS=("PAO1" "PA14" "LESB58")
gene_name="blaNDM-1"
flanking_bp=20000
clinker_distance_cap=50000
threads=16
max_bg_per_country=3
PA445_date="2026"                 # collection date of PA445 (NOT in NCBI metadata)
focal_st=""                       # leave empty to use PA445's ST; or set e.g. "ST308"
beast_clock="strict"             # relaxed | strict
beast_tree_prior="constant"       # constant | skyline
beast_chain_length=100000000
beast_log_every=10000
beast_burnin_pct=10
run_beast="true"                  # set "false" to stop after generating the XML

# Databases
platon_db="${REF_PATH}/platon_db/db/"
plasmidfinder_db="${REF_PATH}/plasmidfinder_db"

# Custom scripts and tools
gbk_region_extraction="${TOOL_PATH}/gbk_region_extraction.py"
prokka_rgi_reconcile="${TOOL_PATH}/prokka_rgi_reconcile.py"
subsample_by_st="${TOOL_PATH}/subsample_by_st.py"
ncbi_accessions_filter="${TOOL_PATH}/ncbi_accessions_filter.py"
rename_contigs="${TOOL_PATH}/rename_contigs.py"
gene_location_filter="${TOOL_PATH}/gene_location_filter.py"
prep_beast_inputs="${TOOL_PATH}/prep_beast_inputs.py"
make_beast_xml="${TOOL_PATH}/make_beast_xml.py"

# Global log
LOG="${COMP_PATH}/PA_comparative.log"
exec > >(tee -a "${LOG}") 2>&1

echo "======================================================"
echo " P. aeruginosa comparative genomics - Started: $(date)"
echo " Focal strain : PA445 (blaNDM-1, chromosomal)"
echo "======================================================"

	############################
	# PHASE 1: DOWNLOAD GENOMES#
	############################

	echo -e "\e[33m ======================================================= \e[0m"
	echo -e "\e[33m PHASE 1: QUERYING AND DOWNLOADING REFERENCE GENOMES      \e[0m"
	echo -e "\e[33m ======================================================= \e[0m"

	# --- 1a. NCBI dataset queries ---
	echo -e "\e[32m Querying NCBI: NDM P. aeruginosa genomes... \e[0m"
	conda run -n ncbi datasets summary genome \
	taxon "Pseudomonas aeruginosa" \
	--assembly-level complete \
	--search "NDM" \
	--limit 30 \
	> "${REFERENCES_PA}/download_logs/summary_NDM.json" \
	2> "${REFERENCES_PA}/download_logs/summary_NDM.stderr.log"

	echo -e "\e[32m Querying NCBI: background P. aeruginosa genomes... \e[0m"
	conda run -n ncbi datasets summary genome \
	taxon "Pseudomonas aeruginosa" \
	--assembly-level complete \
	--limit 500 \
	> "${REFERENCES_PA}/download_logs/summary_background.json"

	# --- 1b. Filter background with ncbi_accessions_filter.py ---
	# NOTE: NDM accessions are ALWAYS retained and do NOT pass through
	# ncbi_accessions_filter.py. The filter script is applied to the
	# background JSON only, then both sets are merged for download.
	echo -e "\e[32m Filtering background genomes (host=human, geography required, dedup-bioproj)... \e[0m"

	conda run -n ncbi python3 "${ncbi_accessions_filter}" \
	"${REFERENCES_PA}/download_logs/summary_NDM.json" \
	"${REFERENCES_PA}/download_logs/summary_background.json" \
	--out-dir "${REFERENCES_PA}/download_logs" \
	--geography \
	--host "Homo sapiens" \
	--dedup-bioproj

	# ncbi_accessions_filter.py is expected to output:
	#   accessions_NDM.txt       - all NDM accessions (passed through unfiltered)
	#   accessions_background.txt - filtered background accessions
	# Merge both into a single master accessions list for bulk download
	ACCESSIONS_NDM="${REFERENCES_PA}/download_logs/accessions_NDM.txt"
	ACCESSIONS_BG="${REFERENCES_PA}/download_logs/accessions_background.txt"
	ACCESSIONS_ALL="${REFERENCES_PA}/download_logs/accessions_all.txt"

	if [[ ! -f "${ACCESSIONS_NDM}" ]]; then
		echo "ERROR: ${ACCESSIONS_NDM} not found - check ncbi_accessions_filter.py output"
		exit 1
	fi
	if [[ ! -f "${ACCESSIONS_BG}" ]]; then
		echo "ERROR: ${ACCESSIONS_BG} not found - check ncbi_accessions_filter.py output"
		exit 1
	fi

	cat "${ACCESSIONS_NDM}" "${ACCESSIONS_BG}" | sort -u > "${ACCESSIONS_ALL}"

	ndm_count=$(wc -l < "${ACCESSIONS_NDM}")
	bg_count=$(wc -l < "${ACCESSIONS_BG}")
	all_count=$(wc -l < "${ACCESSIONS_ALL}")
	echo -e "\e[32m NDM accessions   : ${ndm_count} \e[0m"
	echo -e "\e[32m Background        : ${bg_count} (filtered) \e[0m"
	echo -e "\e[32m Total to download : ${all_count} \e[0m"

	# --- 1c. Build master metadata TSV (NDM + filtered background) ---
	# Parses both JSONs in order; NDM rows labelled "NDM", background "BACKGROUND"
	echo -e "\e[32m Building master metadata TSV... \e[0m"

	python3 << 'PYEOF' > "${REFERENCES_PA}/download_logs/accessions.tsv"

import json, re

ndm_json = "/storage/student9/projects/Paeruginosa/comparative_genomics_Paeruginosa/references/download_logs/summary_NDM.json"
bg_json  = "/storage/student9/projects/Paeruginosa/comparative_genomics_Paeruginosa/references/download_logs/summary_background.json"
bg_acc_file = "/storage/student9/projects/Paeruginosa/comparative_genomics_Paeruginosa/references/download_logs/accessions_background.txt"

with open(ndm_json) as f:
    ndm_data = json.load(f)
with open(bg_json) as f:
    bg_data  = json.load(f)
with open(bg_acc_file) as f:
    bg_keep = set(line.strip() for line in f if line.strip())

seen = set()

def bsm_attr(biosample, key):
    """NCBI datasets v2 nests metadata under biosample['attributes']
    ([{name,value},...]); older/other variants use direct keys. Check both
    so geo_loc_name / host / etc. are not silently reported as 'NA'."""
    if not isinstance(biosample, dict):
        return "NA"
    v = biosample.get(key)
    if v not in (None, "", "NA"):
        return v
    for a in biosample.get("attributes", []) or []:
        if isinstance(a, dict) and str(a.get("name", "")).lower() == key.lower():
            val = a.get("value")
            if val not in (None, "", "NA"):
                return val
    return "NA"

print("\t".join([
    "dataset","accession","strain","assembly_name","assembly_level",
    "contigs","genome_size","contig_n50","geo_location","isolation_source",
    "host","collection_date","submitter","biosample","bioproject"
]))

def process_reports(reports, dataset_name, filter_set=None):
    for r in reports:
        accession = r.get("accession", "NA")
        if accession in seen:
            continue
        # For background: only emit rows whose accession passed the filter
        if filter_set is not None and accession not in filter_set:
            continue
        seen.add(accession)

        strain        = r.get("organism", {}).get("infraspecific_names", {}).get("strain", "NA")
        assembly_info = r.get("assembly_info", {})
        assembly_stats= r.get("assembly_stats", {})
        biosample     = assembly_info.get("biosample", {})

        fields = [
            dataset_name,
            accession,
            strain,
            assembly_info.get("assembly_name", "NA"),
            assembly_info.get("assembly_level", "NA"),
            assembly_stats.get("number_of_component_sequences", "NA"),
            assembly_stats.get("total_sequence_length", "NA"),
            assembly_stats.get("contig_n50", "NA"),
            bsm_attr(biosample, "geo_loc_name"),
            bsm_attr(biosample, "isolation_source"),
            bsm_attr(biosample, "host"),
            bsm_attr(biosample, "collection_date"),
            assembly_info.get("submitter", "NA"),
            biosample.get("accession", "NA"),
            assembly_info.get("bioproject_accession", "NA"),
        ]
        cleaned = [re.sub(r'[\t\n\r]+', ' ', str(x)) if x is not None else "NA" for x in fields]
        print("\t".join(map(str, cleaned)))

# NDM: no filter_set -> all NDM accessions are kept
process_reports(ndm_data.get("reports", []), "NDM", filter_set=None)
# Background: only accessions that passed ncbi_accessions_filter.py
process_reports(bg_data.get("reports",  []), "BACKGROUND", filter_set=bg_keep)

PYEOF

	echo -e "\e[32m Master TSV written: $(tail -n +2 "${REFERENCES_PA}/download_logs/accessions.tsv" | wc -l) genomes \e[0m"

	# --- 1d. Download all selected genomes ---
	echo -e "\e[32m Starting genome downloads... \e[0m"

	tail -n +2 "${REFERENCES_PA}/download_logs/accessions.tsv" | \
while IFS=$'\t' read -r dataset accession strain assembly_name assembly_level \
	contigs genome_size contig_n50 geo_location isolation_source \
	host collection_date submitter biosample bioproject
do
	safe_strain=$(echo "${strain}" | tr ' /\\:*?"<>|' '_')
	if [ "${safe_strain}" = "NA" ] || [ -z "${safe_strain}" ]; then
		safe_strain="${accession}"
	fi

	echo ""
	echo -e "\e[32m Downloading ${safe_strain} [${dataset}] ${accession} \e[0m"

	if [ -f "${REFERENCES_PA}/${safe_strain}.fasta" ]; then
		echo -e "\e[32m   Already exists - skipping \e[0m"
		continue
	fi

	conda run -n ncbi datasets download genome \
	accession "${accession}" \
	--include genome \
	--filename "${REFERENCES_PA}/${safe_strain}.zip" || {
		echo -e "\e[31m   Download failed - skipping \e[0m"
		continue
	}

	if ! python3 -m zipfile -e "${REFERENCES_PA}/${safe_strain}.zip" \
		"${REFERENCES_PA}/${safe_strain}_tmp"; then
		echo -e "\e[31m   Could not unzip ${safe_strain}.zip - skipping \e[0m"
		rm -rf "${REFERENCES_PA}/${safe_strain}_tmp" \
		       "${REFERENCES_PA}/${safe_strain}.zip"
		continue
	fi

	fna_file=$(find "${REFERENCES_PA}/${safe_strain}_tmp" -name "*.fna" | head -1)

	if [ -z "${fna_file}" ]; then
		echo -e "\e[31m   No FASTA found - skipping \e[0m"
		rm -rf "${REFERENCES_PA}/${safe_strain}_tmp" \
		       "${REFERENCES_PA}/${safe_strain}.zip"
		continue
	fi

	cp "${fna_file}" "${REFERENCES_PA}/${safe_strain}.fasta"
	rm -rf "${REFERENCES_PA}/${safe_strain}.zip" \
	       "${REFERENCES_PA}/${safe_strain}_tmp"

	echo -e "\e[32m   Saved: ${safe_strain}.fasta \e[0m"
done

	# --- 1e. Core reference strains + PA445 ---
	echo -e "\e[32m Adding canonical reference strains (PAO1, PA14, LESB58)... \e[0m"

for entry in "GCF_000006765.1:PAO1" "GCF_000014625.1:PA14" "GCF_000026645.1:LESB58"; do
	ref_acc="${entry%%:*}"
	ref_strain="${entry##*:}"

	if [ -f "${REFERENCES_PA}/${ref_strain}.fasta" ]; then
		echo "  ${ref_strain} already exists - skipping"
		continue
	fi

	conda run -n ncbi datasets download genome \
	accession "${ref_acc}" \
	--include genome \
	--filename "${REFERENCES_PA}/${ref_strain}.zip"

	python3 -m zipfile -e "${REFERENCES_PA}/${ref_strain}.zip" \
		"${REFERENCES_PA}/${ref_strain}_tmp"

	find "${REFERENCES_PA}/${ref_strain}_tmp" -name "*.fna" | head -1 | \
		xargs -I{} cp {} "${REFERENCES_PA}/${ref_strain}.fasta"

	rm -rf "${REFERENCES_PA}/${ref_strain}.zip" \
	       "${REFERENCES_PA}/${ref_strain}_tmp"
	echo -e "\e[32m   Saved: ${ref_strain}.fasta \e[0m"
done

	cp "${PA445}" "${REFERENCES_PA}/PA445.fasta"

	echo ""
	echo -e "\e[32m All genomes ready ($(ls "${REFERENCES_PA}/"*.fasta | wc -l) total): \e[0m"
	ls -lh "${REFERENCES_PA}/"*.fasta

	#############################################
	# PHASE 1 ANALYSIS: MLST / PASTY / RGI      #
	#                   PROKKA / PANAROO        #
	# Run on FULL dataset before subsampling    #
	#############################################

	echo -e "\e[33m ===================================================== \e[0m"
	echo -e "\e[33m PHASE 1 ANALYSIS: FULL DATASET - MLST / PASTY / RGI   \e[0m"
	echo -e "\e[33m                   PROKKA / PANAROO                    \e[0m"
	echo -e "\e[33m ===================================================== \e[0m"

	# -- MLST --------------------------------------------------------------

	echo -e "\e[31m ================= \e[0m"
	echo -e "\e[31m MLST: ALL STRAINS \e[0m"
	echo -e "\e[31m ================= \e[0m"

	conda run -n BPtyping mlst \
	--scheme paeruginosa \
	--threads "${threads}" \
	"${REFERENCES_PA}/"*.fasta \
	> "${MLST_PA}/all.mlst.tsv"

	echo -e "\e[32m MLST results: \e[0m"
	cat "${MLST_PA}/all.mlst.tsv"

	# -- PASTY -------------------------------------------------------------

	echo -e "\e[31m ================== \e[0m"
	echo -e "\e[31m PASTY: ALL STRAINS \e[0m"
	echo -e "\e[31m ================== \e[0m"

for fasta in "${REFERENCES_PA}/"*.fasta;
do
	strain=$(basename "${fasta}" .fasta)

	conda run -n BPtyping pasty \
	--force \
	--prefix "${strain}.pasty" \
	--input "${fasta}" \
	--outdir "${PASTY_PA}"
done

	# Merge pasty outputs
	head -n 1 "${PASTY_PA}/$(ls "${PASTY_PA}/"*.pasty.tsv | head -1 | xargs basename)" \
		> "${PASTY_PA}/PA_serotypes.tsv"
for fasta in "${REFERENCES_PA}/"*.fasta;
do
	strain=$(basename "${fasta}" .fasta)

	if [ -f "${PASTY_PA}/${strain}.pasty.tsv" ]; then
		tail -n +2 "${PASTY_PA}/${strain}.pasty.tsv" >> "${PASTY_PA}/PA_serotypes.tsv"
	fi
done
	echo -e "\e[32m Pasty serotypes: ${PASTY_PA}/PA_serotypes.tsv \e[0m"

	# -- RGI (full dataset) ------------------------------------------------

	echo -e "\e[31m ================= \e[0m"
	echo -e "\e[31m RGI: FULL DATASET         \e[0m"
	echo -e "\e[31m ================= \e[0m"

	echo -e "\e[32m Cleaning FASTA headers for RGI compatibility... \e[0m"

	echo -e "\e[32m Renaming contigs to RGI/Prodigal-safe names (bare contigNNNNN + one shared map)... \e[0m"

# One combined map for the whole full-dataset run: sample -> contig -> original header.
RGI_PA_MAP="${RGI_PA}/rgi_inputs/contig_name_map.tsv"
: > "${RGI_PA_MAP}"

for fasta in "${REFERENCES_PA}/"*.fasta;
do
	strain=$(basename "${fasta}" .fasta)

	output="${RGI_PA}/rgi_inputs/${strain}.fasta"
	# Contigs -> c00001, c00002, ... (no underscore; one sequence line each so the
	# faidx 'inconsistent line length' failure - the real cause of RGI's
	# "Requested rname ... does not exist" - cannot happen). --validate builds the
	# .fai with htslib (what RGI uses) and fails loudly if the file is unindexable.
	conda run -n rgi_env python3 "${rename_contigs}" \
	"${fasta}" "${output}" \
	--sample "${strain}" \
	--map "${RGI_PA_MAP}" \
	--append \
	--prefix c \
	--width 0 \
	--validate
done

for fasta in "${RGI_PA}/rgi_inputs/"*.fasta;
do
	strain=$(basename "${fasta}" .fasta)

	echo -e "\e[32m Running RGI on ${strain}... \e[0m"
	conda run -n rgi_env rgi main \
	--input_sequence "${fasta}" \
	--output_file "${RGI_PA}/${strain}_rgi" \
	--input_type contig \
	--alignment_tool DIAMOND \
	--num_threads "${threads}" \
	--clean
done

	echo -e "\e[32m Generating RGI heatmap for subsampled dataset... \e[0m"
	conda run -n rgi_env rgi heatmap \
	--input "${RGI_PA}/" \
	--output "${RGI_PA}/PA445.AMR_heatmap_full"

	echo -e "\e[32m RGI full-dataset complete: ${RGI_PA}/ \e[0m"

	# -- PROKKA (full dataset) ---------------------------------------------

	echo -e "\e[31m ==================================== \e[0m"
	echo -e "\e[31m PROKKA: ANNOTATING ALL STRAINS       \e[0m"
	echo -e "\e[31m ==================================== \e[0m"

	conda run -n BPannotation prokka \
	--force \
	--cpus "${threads}" \
	--genus Pseudomonas \
	--species aeruginosa \
	--prefix PA445.prokka \
	--outdir "${PROKKA_PA}/PA445.prokka" \
	"${PA445}"

	cp "${PROKKA_PA}/PA445.prokka/PA445.prokka.gff" \
	   "${PANAROO_PA}/gff_inputs/PA445.gff"

for fasta in "${REFERENCES_PA}/"*.fasta;
do
	strain=$(basename "${fasta}" .fasta)
	[ "${strain}" == "PA445" ] && continue

	mkdir -p "${PROKKA_PA}/${strain}.prokka"
	echo -e "\e[32m Annotating ${strain} with Prokka... \e[0m"

	conda run -n BPannotation prokka \
	--force \
	--cpus "${threads}" \
	--genus Pseudomonas \
	--species aeruginosa \
	--prefix "${strain}.prokka" \
	--outdir "${PROKKA_PA}/${strain}.prokka" \
	"${fasta}"

	cp "${PROKKA_PA}/${strain}.prokka/${strain}.prokka.gff" \
	   "${PANAROO_PA}/gff_inputs/${strain}.gff"
	echo -e "\e[32m   ${strain} annotation complete \e[0m"
done

	# -- PANAROO (full dataset) --------------------------------------------

	echo -e "\e[31m ================================ \e[0m"
	echo -e "\e[31m PANAROO: FULL DATASET PAN-GENOME \e[0m"
	echo -e "\e[31m ================================ \e[0m"

	echo -e "\e[32m GFF files for pangenome: $(ls "${PANAROO_PA}/gff_inputs/"*.gff | wc -l) \e[0m"

	conda run -n pangenome panaroo \
	--threads "${threads}" \
	--input "${PANAROO_PA}/gff_inputs/"*.gff \
	--out_dir "${PANAROO_PA}/results" \
	--clean-mode moderate \
	--alignment core \
	--core_threshold 0.90 \
	--aligner mafft

	echo -e "\e[32m Pan-genome complete: ${PANAROO_PA}/results/ \e[0m"
	if [ -f "${PANAROO_PA}/results/summary_statistics.txt" ]; then
		echo -e "\e[32m Pan-genome summary: \e[0m"
		cat "${PANAROO_PA}/results/summary_statistics.txt"
	fi

	# -- FASTANI (full dataset) ---------------------------------------------

	echo -e "\e[31m ===================== \e[0m"
	echo -e "\e[31m FASTANI: FULL DATASET \e[0m"
	echo -e "\e[31m ===================== \e[0m"

	ls "${REFERENCES_PA}/"*.fasta > "${ANI_PA}/full/genome_list.txt"

	conda run -n ani fastANI \
	--ql "${ANI_PA}/full/genome_list.txt" \
	--rl "${ANI_PA}/full/genome_list.txt" \
	--threads "${threads}" \
	-o "${ANI_PA}/full/PA445.fastANI_full.tsv"

	echo -e "\e[32m fastANI (full dataset): ${ANI_PA}/full/PA445.fastANI_full.tsv \e[0m"

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m ANICLUSTERMAP: FULL DATASET \e[0m"
	echo -e "\e[31m =========================== \e[0m"

	conda run -n ani ANIclustermap \
	-i "${REFERENCES_PA}" \
	-o "${ANI_PA}/full" \
	--fig_width 15 \
	--fig_height 12 

	echo -e "\e[32m ANI heatmap (full dataset): ${ANI_PA}/full/ANIclustermap.pdf \e[0m"

	####################################################
	# PHASE 2: IDENTIFY TRUE NDM-1 STRAINS             #
	# prokka_rgi_reconcile.py across the full dataset  #
	####################################################

	echo -e "\e[33m ========================================================= \e[0m"
	echo -e "\e[33m PHASE 2: PROKKA + RGI RECONCILIATION - TRUE NDM-1 STRAINS \e[0m"
	echo -e "\e[33m ========================================================= \e[0m"

	rm -rf "${CLINKER_PA}/truepos_gbk"
	mkdir -p "${CLINKER_PA}/truepos_gbk"

	conda run -n phylogeny python3 "${prokka_rgi_reconcile}" \
	--prokka-gene "blaNDM-1" \
	--prokka-inputs "${PROKKA_PA}" \
	--rgi-gene "NDM-1" \
	--rgi-inputs "${RGI_PA}" \
	--rgi-cutoff Strict,Perfect \
	--truepos-prokka-paths \
	--outdir "${CLINKER_PA}/truepos_gbk"

	TRUEPOS_PATHS="${CLINKER_PA}/truepos_gbk/blaNDM1.truepos_prokka_gbk_paths.txt"

	if [[ ! -f "${TRUEPOS_PATHS}" ]]; then
		echo "ERROR: True-positive paths file not found: ${TRUEPOS_PATHS}"
		exit 1
	fi

	truepos_count=$(wc -l < "${TRUEPOS_PATHS}")
	echo -e "\e[32m True NDM-1 strains confirmed: ${truepos_count} \e[0m"
	cat "${TRUEPOS_PATHS}"

	# Build a plain list of true-NDM-1 strain names (for downstream use)
	NDM1_STRAINS="${CLINKER_PA}/truepos_gbk/ndm1_strain_names.txt"
	while IFS= read -r gbk_path; do
		[[ -z "${gbk_path}" ]] && continue
		basename "${gbk_path}" .prokka.gbk
	done < "${TRUEPOS_PATHS}" > "${NDM1_STRAINS}"
	echo -e "\e[32m NDM-1 strain name list: ${NDM1_STRAINS} \e[0m"

	#######################################################
	# PHASE 3: SUBSAMPLING (ST-aware)                     #
	#   - All confirmed NDM-1 strains (protected)         #
	#   - 1 background rep per ST (clonal dedup),         #
	#     a few for STs sharing an NDM-1 lineage          #
	#   - then <=max_bg_per_country per country, dropping #
	#     ST-redundant genomes first                      #
	#######################################################

	echo -e "\e[33m ==================================================================== \e[0m"
	echo -e "\e[33m PHASE 3: SUBSAMPLING DATASET (ST-aware)                              \e[0m"
	echo -e "\e[33m All NDM-1 strains + 1/ST background, <=${max_bg_per_country}/country \e[0m"
	echo -e "\e[33m ==================================================================== \e[0m"

	SUBSAMPLE_LIST="${SUBSAMPLE_PA}/subsampled_strains.txt"

	# Dedup background by sequence type (removes clonal over-representation, which
	# country caps alone cannot), keep the best assembly per ST, give STs that share
	# a lineage with the NDM-1 strains a couple of context reps, then apply an
	# ST-aware per-country cap. PA445 is protected via the NDM-1 list. References
	# (PAO1/PA14/LESB58) are added later in the Phase 7 controlled view, not here.
	conda run -n phylogeny python3 "${subsample_by_st}" \
	--metadata "${REFERENCES_PA}/download_logs/accessions.tsv" \
	--mlst "${MLST_PA}/all.mlst.tsv" \
	--protect "${NDM1_STRAINS}" \
	--max-per-st 1 \
	--context-per-st 2 \
	--max-per-country "${max_bg_per_country}" \
	--out-strains "${SUBSAMPLE_LIST}" \
	--out-table "${SUBSAMPLE_PA}/subsample_audit.tsv"

	sub_count=$(wc -l < "${SUBSAMPLE_LIST}")
	echo -e "\e[32m Subsampled strain list written: ${sub_count} strains -> ${SUBSAMPLE_LIST} \e[0m"

	# Copy subsampled FASTAs to dedicated directory
	echo -e "\e[32m Copying subsampled FASTAs... \e[0m"
	while IFS= read -r strain_name; do
		[[ -z "${strain_name}" ]] && continue
		src="${REFERENCES_PA}/${strain_name}.fasta"
		if [ -f "${src}" ]; then
			cp "${src}" "${SUBSAMPLE_FASTA}/${strain_name}.fasta"
		else
			echo -e "\e[31m   WARNING: ${src} not found - skipping \e[0m"
		fi
	done < "${SUBSAMPLE_LIST}"

	echo -e "\e[32m Subsampled FASTAs: $(ls "${SUBSAMPLE_FASTA}/"*.fasta 2>/dev/null | wc -l) \e[0m"
	ls -lh "${SUBSAMPLE_FASTA}/"

	######################################################
	# PHASE 4a: RGI HEATMAP - SUBSAMPLED DATASET         #
	# Reuse existing full-dataset RGI JSON outputs where #
	# available; only re-run RGI if output is missing.   #
	######################################################

	echo -e "\e[33m ================================================== \e[0m"
	echo -e "\e[33m PHASE 4a: RGI HEATMAP - SUBSAMPLED DATASET         \e[0m"
	echo -e "\e[33m (reusing full-dataset RGI outputs where available) \e[0m"
	echo -e "\e[33m ================================================== \e[0m"

# Shared map for any contigs freshly renamed in the subsampled re-run
# (most strains reuse the full-dataset RGI output and its map instead).
RGI_SUB_MAP="${RGI_SUB}/rgi_inputs/contig_name_map.tsv"
: > "${RGI_SUB_MAP}"

while IFS= read -r strain_name; do
	[[ -z "${strain_name}" ]] && continue

	src_txt="${RGI_PA}/${strain_name}_rgi.txt"
	src_json="${RGI_PA}/${strain_name}_rgi.json"
	dst_txt="${RGI_SUB}/${strain_name}_rgi.txt"
	dst_json="${RGI_SUB}/${strain_name}_rgi.json"

	if [ -f "${src_txt}" ] && [ -f "${src_json}" ]; then
		# Reuse existing RGI outputs (symlink avoids disk duplication)
		ln -sf "${src_txt}"  "${dst_txt}"
		ln -sf "${src_json}" "${dst_json}"
		echo -e "\e[32m   ${strain_name}: reusing existing RGI output \e[0m"
	else
		# Run RGI if missing (e.g. core references or new strains)
		echo -e "\e[32m   ${strain_name}: running RGI (output not found in full-dataset dir) \e[0m"

		fasta="${RGI_SUB}/rgi_inputs/${strain_name}.fasta"
		if [ ! -f "${fasta}" ]; then
			# Rename contigs (same RGI/Prodigal-safe scheme as the full run)
			src_fasta="${REFERENCES_PA}/${strain_name}.fasta"
			if [ -f "${src_fasta}" ]; then
				conda run -n rgi_env python3 "${rename_contigs}" \
				"${src_fasta}" "${fasta}" \
				--sample "${strain_name}" \
				--map "${RGI_SUB_MAP}" \
				--append \
				--prefix c \
				--width 0 \
				--validate
			else
				echo -e "\e[31m   WARNING: source FASTA not found for ${strain_name} - skipping RGI \e[0m"
				continue
			fi
		fi

		conda run -n rgi_env rgi main \
		--input_sequence "${fasta}" \
		--output_file "${RGI_SUB}/${strain_name}_rgi" \
		--input_type contig \
		--alignment_tool DIAMOND \
		--num_threads "${threads}" \
		--clean
	fi
done < "${SUBSAMPLE_LIST}"

	echo -e "\e[32m Generating RGI heatmap for subsampled dataset... \e[0m"
	conda run -n rgi_env rgi heatmap \
	--input "${RGI_SUB}/" \
	--output "${RGI_SUB}/PA445.AMR_heatmap_subsampled"

	echo -e "\e[32m AMR heatmap (subsampled): ${RGI_SUB}/PA445.AMR_heatmap_subsampled \e[0m"

	######################################################
	# PHASE 4b: IQ-TREE PHYLOGENY - SUBSAMPLED DATASET   #
	######################################################

	echo -e "\e[33m ================================================ \e[0m"
	echo -e "\e[33m PHASE 4b: IQ-TREE PHYLOGENY - SUBSAMPLED DATASET \e[0m"
	echo -e "\e[33m ================================================ \e[0m"

	# Build subsampled GFF list and run a targeted Panaroo for core alignment
	PANAROO_SUB="${COMP_PATH}/panaroo_subsampled"
	mkdir -p "${PANAROO_SUB}/gff_inputs"
	mkdir -p "${PANAROO_SUB}/results"

while IFS= read -r strain_name; do
	[[ -z "${strain_name}" ]] && continue
	src_gff="${PANAROO_PA}/gff_inputs/${strain_name}.gff"
	if [ -f "${src_gff}" ]; then
		ln -sf "${src_gff}" "${PANAROO_SUB}/gff_inputs/${strain_name}.gff"
	else
		echo -e "\e[31m   WARNING: GFF not found for ${strain_name} - will be excluded from phylogeny \e[0m"
	fi
done < "${SUBSAMPLE_LIST}"

	gff_sub_count=$(ls "${PANAROO_SUB}/gff_inputs/"*.gff 2>/dev/null | wc -l)
	echo -e "\e[32m GFF files for subsampled panaroo: ${gff_sub_count} \e[0m"

	conda run -n pangenome panaroo \
	--threads "${threads}" \
	--input "${PANAROO_SUB}/gff_inputs/"*.gff \
	--out_dir "${PANAROO_SUB}/results" \
	--clean-mode moderate \
	--alignment core \
	--core_threshold 0.90 \
	--aligner mafft

	if [ ! -f "${PANAROO_SUB}/results/core_gene_alignment_filtered.aln" ]; then
		echo "ERROR: subsampled core alignment not found - check panaroo_subsampled/results/"
		exit 1
	fi

	conda run -n phylogeny iqtree \
	-s "${PANAROO_SUB}/results/core_gene_alignment_filtered.aln" \
	-m GTR+F+G4 \
	-bb 1000 \
	-nt "${threads}" \
	--prefix "${IQTREE_SUB}/PA445.phylogeny_subsampled"

	echo -e "\e[32m Phylogeny (subsampled): ${IQTREE_SUB}/PA445.phylogeny_subsampled.treefile \e[0m"

	# -- FASTANI (subsampled dataset) ----------------------------------------

	echo -e "\e[31m ============================ \e[0m"
	echo -e "\e[31m FASTANI: SUBSAMPLED DATASET \e[0m"
	echo -e "\e[31m ============================ \e[0m"

	ls "${SUBSAMPLE_FASTA}/"*.fasta > "${ANI_PA}/subsampled/genome_list.txt"

	conda run -n ani fastANI \
	--ql "${ANI_PA}/subsampled/genome_list.txt" \
	--rl "${ANI_PA}/subsampled/genome_list.txt" \
	--threads "${threads}" \
	-o "${ANI_PA}/subsampled/PA445.fastANI_subsampled.tsv"

	echo -e "\e[32m fastANI (subsampled dataset): ${ANI_PA}/subsampled/PA445.fastANI_subsampled.tsv \e[0m"

	echo -e "\e[31m ================================= \e[0m"
	echo -e "\e[31m ANICLUSTERMAP: SUBSAMPLED DATASET \e[0m"
	echo -e "\e[31m ================================= \e[0m"

	conda run -n ani ANIclustermap \
	-i "${SUBSAMPLE_FASTA}" \
	-o "${ANI_PA}/subsampled" \
	--fig_width 15 \
	--fig_height 12 \

	echo -e "\e[32m ANI heatmap (subsampled dataset): ${ANI_PA}/subsampled/aniclustermap/ANIclustermap.pdf \e[0m"

	########################################################
	# PHASE 5: GBK REGION EXTRACTION + CLINKER             #
	# NDM-1 confirmed strains only                         #
	########################################################

	echo -e "\e[33m ===================================================== \e[0m"
	echo -e "\e[33m PHASE 5: GBK EXTRACTION + CLINKER (NDM-1 STRAINS ONLY) \e[0m"
	echo -e "\e[33m ===================================================== \e[0m"

	rm -rf "${CLINKER_PA}/gbk_inputs"
	mkdir -p "${CLINKER_PA}/gbk_inputs"

while IFS= read -r full_gbk; do
	[[ -z "${full_gbk}" ]] && continue
	[[ -f "${full_gbk}" ]] || { echo "WARNING: GBK not found: ${full_gbk} - skipping"; continue; }

	strain=$(basename "${full_gbk}" .prokka.gbk)
	echo -e "\e[32m Extracting ${gene_name} region for ${strain}... \e[0m"

	conda run -n phylogeny python3 "${gbk_region_extraction}" \
	"${full_gbk}" \
	"${CLINKER_PA}/gbk_inputs" \
	"${strain}" \
	--gene "${gene_name}" \
	--flank "${flanking_bp}" \
	${clinker_distance_cap:+--distance-cap "${clinker_distance_cap}"}

done < "${TRUEPOS_PATHS}"

	echo -e "\e[32m Region GBK extraction complete: \e[0m"
	ls -lh "${CLINKER_PA}/gbk_inputs/"

	gbk_count=$(ls "${CLINKER_PA}/gbk_inputs/"*.gbk 2>/dev/null | wc -l)
	echo -e "\e[32m GBK files available for synteny: ${gbk_count} \e[0m"

	if [ "${gbk_count}" -gt 0 ]; then
		conda run -n phylogeny clinker \
		"${CLINKER_PA}/gbk_inputs/"*.gbk \
		--plot "${CLINKER_PA}/blaNDM1.clinker.html" \
		--identity 0.3

		echo -e "\e[32m Synteny plot: ${CLINKER_PA}/blaNDM1.clinker.html \e[0m"
	else
		echo -e "\e[33m WARNING: No GBK files found - skipping clinker \e[0m"
	fi

	######################################################
	# PHASE 6: PLASMID ANALYSIS - NDM-1 STRAINS ONLY     #
	######################################################

	echo -e "\e[33m ===================================================== \e[0m"
	echo -e "\e[33m PHASE 6: PLASMID ANALYSIS (NDM-1 STRAINS ONLY)        \e[0m"
	echo -e "\e[33m   Rationale: supervisor priority is chromosomal       \e[0m"
	echo -e "\e[33m   location of blaNDM-1, as seen in PA445              \e[0m"
	echo -e "\e[33m ===================================================== \e[0m"

while IFS= read -r strain_name; do
	[[ -z "${strain_name}" ]] && continue

	fasta="${REFERENCES_PA}/${strain_name}.fasta"
	[ "${strain_name}" == "PA445" ] && fasta="${PA445}"

	if [ ! -f "${fasta}" ]; then
		echo -e "\e[31m   WARNING: FASTA not found for ${strain_name} - skipping plasmid analysis \e[0m"
		continue
	fi

	mkdir -p "${PLATON_PA}/${strain_name}.platon"
	mkdir -p "${PLASMIDFINDER_PA}/${strain_name}.plasmidfinder"
	mkdir -p "${MOBSUITE_PA}/${strain_name}.mobsuite"

	PLATON_OUT="${PLATON_PA}/${strain_name}.platon"
	PLASMIDFINDER_OUT="${PLASMIDFINDER_PA}/${strain_name}.plasmidfinder"
	MOBSUITE_OUT="${MOBSUITE_PA}/${strain_name}.mobsuite"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m PLATON: ${strain_name} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	conda run -n plasmid platon \
	--db "${platon_db}" \
	--output "${PLATON_OUT}" \
	--prefix "${strain_name}" \
	--mode sensitivity \
	--threads "${threads}" \
	"${fasta}"

	echo -e "\e[32m   Platon complete for ${strain_name} \e[0m"

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m PLASMIDFINDER: ${strain_name} \e[0m"
	echo -e "\e[31m =========================== \e[0m"

	conda run -n plasmid plasmidfinder.py \
	-i "${fasta}" \
	-o "${PLASMIDFINDER_OUT}" \
	-p "${plasmidfinder_db}" \
	-l 0.60 \
	-t 0.80 \
	-x

	echo -e "\e[32m   PlasmidFinder complete for ${strain_name} \e[0m"

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m MOB-SUITE: ${strain_name} \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	conda run -n plasmid mob_recon \
	--infile "${fasta}" \
	--outdir "${MOBSUITE_OUT}" \
	--num_threads "${threads}" \
	--force

	echo -e "\e[32m   MOB-suite complete for ${strain_name} \e[0m"
	echo -e "\e[32m     chromosome.fasta    - chromosomal contigs \e[0m"
	echo -e "\e[32m     plasmid_*.fasta     - reconstructed plasmids \e[0m"
	echo -e "\e[32m     mobtyper_results.txt - plasmid mobility/typing \e[0m"

done < "${NDM1_STRAINS}"

	#############################################################
	# PHASE 7: CHROMOSOMAL-NDM-1 CONTROLLED VIEW                #
	#   Reuse RGI + redraw phylogeny on a focused set:          #
	#     - background : <=5 per country (same as Phase 3)      #
	#     - references : PAO1, PA14, LESB58 (phylogeny anchors) #
	#     - NDM-1 strains with blaNDM-1 ON THE CHROMOSOME       #
	#     - focal isolate PA445                                 #
	#   (Phases 1-6 above are untouched; this block is additive)#
	#############################################################

	echo -e "\e[33m ========================================== \e[0m"
	echo -e "\e[33m PHASE 7: CHROMOSOMAL-NDM-1 CONTROLLED VIEW \e[0m"
	echo -e "\e[33m ========================================== \e[0m"

	# --- 7a. Which NDM-1 strains carry blaNDM-1 on the chromosome? ----------
	# Location is read from the plasmid tools (MOB-suite/Platon, run in Phase 6),
	# not re-derived from RGI: the NDM-1 calls are already true-positive-verified
	# by prokka_rgi_reconcile.py, so only chromosome-vs-plasmid remains. Prokka
	# tells us which contig carries the gene; MOB-suite/Platon classify that contig.
	echo -e "\e[32m Classifying ${gene_name} location (chromosome vs plasmid) from Prokka + MOB-suite/Platon... \e[0m"

	CHROM_NDM_STRAINS="${SUBSUB_PA}/chromosomal_ndm1_strains.txt"
	NDM1_LOCATION_REPORT="${SUBSUB_PA}/ndm1_location_report.tsv"
	CHROM_NDM_FASTAS="${SUBSUB_PA}/chromosomal_ndm1_fastas.txt"

	conda run -n phylogeny python3 "${gene_location_filter}" \
	--gene "${gene_name}" \
	--prokka-dir "${PROKKA_PA}" \
	--strains "${NDM1_STRAINS}" \
	--mobsuite-dir "${MOBSUITE_PA}" \
	--platon-dir "${PLATON_PA}" \
	--location-source consensus \
	--keep chromosome \
	--fasta-dir "${REFERENCES_PA}" \
	--out-table "${NDM1_LOCATION_REPORT}" \
	--out-strains "${CHROM_NDM_STRAINS}" \
	--out-fasta-list "${CHROM_NDM_FASTAS}"

	chrom_ndm_count=$(wc -l < "${CHROM_NDM_STRAINS}")
	echo -e "\e[32m Chromosomal ${gene_name} strains: ${chrom_ndm_count} (audit: ${NDM1_LOCATION_REPORT}) \e[0m"

	# --- 7a-2. Strains sharing PA445's exact ST -----------------------------
	# Supervisor request: only ~5 global isolates (out of the top-500 NCBI
	# screen) share PA445's ST. These are kept regardless of country caps,
	# since they're the closest global relatives available for this ST.
	echo -e "\e[31m ======================================== \e[0m"
	echo -e "\e[31m Identifying strains sharing PA445's ST   \e[0m"
	echo -e "\e[31m ======================================== \e[0m"

	SAME_ST_STRAINS="${SUBSUB_PA}/same_st_as_PA445.txt"

	# FILE column holds full paths (e.g. .../PA445.fasta) - strip to basename
	PA445_ST=$(awk -F'\t' '{n=split($1,a,"/"); f=a[n]; sub(/\.fasta$/,"",f); if (f=="PA445") print $3}' \
		"${MLST_PA}/all.mlst.tsv")

	if [[ -z "${PA445_ST}" || "${PA445_ST}" == "-" ]]; then
		echo -e "\e[32m WARNING: Could not determine PA445 ST from ${MLST_PA}/all.mlst.tsv - skipping ST-mate merge \e[0m"
		: > "${SAME_ST_STRAINS}"
	else
		echo -e "\e[32m PA445 ST: ${PA445_ST} \e[0m"
		awk -F'\t' -v st="${PA445_ST}" \
		'{n=split($1,a,"/"); f=a[n]; sub(/\.fasta$/,"",f); if ($3==st && f!="PA445") print f}' \
		"${MLST_PA}/all.mlst.tsv" > "${SAME_ST_STRAINS}"

		echo -e "\e[32m Strains sharing PA445's ST: $(wc -l < "${SAME_ST_STRAINS}") \e[0m"
		cat "${SAME_ST_STRAINS}"
	fi

	# --- 7b. Assemble the controlled strain set -----------------------------
	# Background = Phase-3 subsample minus NDM-1 strains => identical <=3/country set.
	SUBSUB_LIST="${SUBSUB_PA}/subsubsampled_strains.txt"
	BG_ONLY="${SUBSUB_PA}/background_only.txt"
	comm -23 <(sort -u "${SUBSAMPLE_LIST}") <(sort -u "${NDM1_STRAINS}") > "${BG_ONLY}"

	{
		cat "${BG_ONLY}"
		printf '%s\n' "${REF_ANCHORS[@]}"
		cat "${CHROM_NDM_STRAINS}"
		cat "${SAME_ST_STRAINS}"
		echo "PA445"
	} | sort -u > "${SUBSUB_LIST}"

	subsub_count=$(wc -l < "${SUBSUB_LIST}")
	echo -e "\e[32m Controlled set: ${subsub_count} strains -> ${SUBSUB_LIST} \e[0m"
	echo -e "\e[32m   background: $(wc -l < "${BG_ONLY}"), references: ${#REF_ANCHORS[@]}, chromosomal NDM-1: ${chrom_ndm_count}, ST-mates: $(wc -l < "${SAME_ST_STRAINS}"), +PA445 \e[0m"

	# Copy controlled-set FASTAs to a dedicated directory (PA445 from its source)
	while IFS= read -r strain_name; do
		[[ -z "${strain_name}" ]] && continue
		if [ "${strain_name}" = "PA445" ]; then
			src="${PA445}"
		else
			src="${REFERENCES_PA}/${strain_name}.fasta"
		fi
		if [ -f "${src}" ]; then
			cp "${src}" "${SUBSUB_FASTA}/${strain_name}.fasta"
		else
			echo -e "\e[31m   WARNING: ${src} not found - skipping copy \e[0m"
		fi
	done < "${SUBSUB_LIST}"

	# --- 7c. RGI: reuse full-dataset outputs, heatmap -----------------------
	echo -e "\e[31m ================================== \e[0m"
	echo -e "\e[31m RGI heatmap: Subsubsampled dataset \e[0m"
	echo -e "\e[31m ================================== \e[0m"

RGI_SUBSUB_MAP="${RGI_SUBSUB}/rgi_inputs/contig_name_map.tsv"
: > "${RGI_SUBSUB_MAP}"

while IFS= read -r strain_name; do
	[[ -z "${strain_name}" ]] && continue
	src_txt="${RGI_PA}/${strain_name}_rgi.txt"
	src_json="${RGI_PA}/${strain_name}_rgi.json"
	dst_txt="${RGI_SUBSUB}/${strain_name}_rgi.txt"
	dst_json="${RGI_SUBSUB}/${strain_name}_rgi.json"

	if [ -f "${src_txt}" ] && [ -f "${src_json}" ]; then
		ln -sf "${src_txt}"  "${dst_txt}"
		ln -sf "${src_json}" "${dst_json}"
		echo -e "\e[32m   ${strain_name}: reusing existing RGI output \e[0m"
	else
		echo -e "\e[32m   ${strain_name}: RGI output missing in full set - running RGI \e[0m"
		fasta="${RGI_SUBSUB}/rgi_inputs/${strain_name}.fasta"
		if [ "${strain_name}" = "PA445" ]; then
			src_fasta="${PA445}"
		else
			src_fasta="${REFERENCES_PA}/${strain_name}.fasta"
		fi
		if [ -f "${src_fasta}" ]; then
			conda run -n rgi_env python3 "${rename_contigs}" \
			"${src_fasta}" "${fasta}" \
			--sample "${strain_name}" \
			--map "${RGI_SUBSUB_MAP}" \
			--append \
			--prefix c \
			--width 0 \
			--validate

			conda run -n rgi_env rgi main \
			--input_sequence "${fasta}" \
			--output_file "${RGI_SUBSUB}/${strain_name}_rgi" \
			--input_type contig \
			--alignment_tool DIAMOND \
			--num_threads "${threads}" \
			--clean
		else
			echo -e "\e[32m   WARNING: source FASTA not found for ${strain_name} - skipping RGI \e[0m"
			continue
		fi
	fi
done < "${SUBSUB_LIST}"

	echo -e "\e[32m Generating RGI heatmap for controlled set... \e[0m"
	conda run -n rgi_env rgi heatmap \
	--input "${RGI_SUBSUB}/" \
	--output "${RGI_SUBSUB}/PA445.AMR_heatmap_chromosomal_ndm1"

	echo -e "\e[32m AMR heatmap (controlled): ${RGI_SUBSUB}/PA445.AMR_heatmap_chromosomal_ndm1 \e[0m"

	# --- 7d. Phylogeny on the controlled set --------------------------------
	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m IQ-TREE: Controlled set \e[0m"
	echo -e "\e[31m ======================= \e[0m"

while IFS= read -r strain_name; do
	[[ -z "${strain_name}" ]] && continue
	src_gff="${PANAROO_PA}/gff_inputs/${strain_name}.gff"
	if [ -f "${src_gff}" ]; then
		ln -sf "${src_gff}" "${PANAROO_SUBSUB}/gff_inputs/${strain_name}.gff"
	else
		echo -e "\e[31m   WARNING: GFF not found for ${strain_name} - excluded from phylogeny \e[0m"
	fi
done < "${SUBSUB_LIST}"

	gff_subsub_count=$(ls "${PANAROO_SUBSUB}/gff_inputs/"*.gff 2>/dev/null | wc -l)
	echo -e "\e[32m GFF files for controlled panaroo: ${gff_subsub_count} \e[0m"

	conda run -n pangenome panaroo \
	--threads "${threads}" \
	--input "${PANAROO_SUBSUB}/gff_inputs/"*.gff \
	--out_dir "${PANAROO_SUBSUB}/results" \
	--clean-mode sensitive \
	--alignment core \
	--core_threshold 0.90 \
	--aligner mafft

	if [ ! -f "${PANAROO_SUBSUB}/results/core_gene_alignment_filtered.aln" ]; then
		echo "ERROR: controlled-set core alignment not found - check panaroo_subsubsampled/results/"
		exit 1
	fi

	conda run -n phylogeny iqtree \
	-s "${PANAROO_SUBSUB}/results/core_gene_alignment_filtered.aln" \
	-m GTR+F+G4 \
	-bb 1000 \
	-nt "${threads}" \
	--prefix "${IQTREE_SUBSUB}/PA445.phylogeny_chromosomal_ndm1"

	echo -e "\e[32m Phylogeny (controlled): ${IQTREE_SUBSUB}/PA445.phylogeny_chromosomal_ndm1.treefile \e[0m"

	# -- FASTANI (subsubsampled/controlled dataset) --------------------------

	echo -e "\e[31m =============================== \e[0m"
	echo -e "\e[31m FASTANI: SUBSUBSAMPLED DATASET \e[0m"
	echo -e "\e[31m =============================== \e[0m"

	ls "${SUBSUB_FASTA}/"*.fasta > "${ANI_PA}/subsubsampled/genome_list.txt"

	conda run -n ani fastANI \
	--ql "${ANI_PA}/subsubsampled/genome_list.txt" \
	--rl "${ANI_PA}/subsubsampled/genome_list.txt" \
	--threads "${threads}" \
	-o "${ANI_PA}/subsubsampled/PA445.fastANI_subsubsampled.tsv"

	echo -e "\e[32m fastANI (subsubsampled dataset): ${ANI_PA}/subsubsampled/PA445.fastANI_subsubsampled.tsv \e[0m"

	echo -e "\e[31m ==================================== \e[0m"
	echo -e "\e[31m ANICLUSTERMAP: SUBSUBSAMPLED DATASET \e[0m"
	echo -e "\e[31m ==================================== \e[0m"

	conda run -n ani ANIclustermap \
	-i "${SUBSUB_FASTA}" \
	-o "${ANI_PA}/subsubsampled" \
	--fig_width 12 \
	--fig_height 10 

	echo -e "\e[32m ANI heatmap (subsubsampled dataset): ${ANI_PA}/subsubsampled/aniclustermap/ANIclustermap.pdf \e[0m"

	#############################################################
	# PHASE 8: TIME-CALIBRATED PHYLOGENY (BEAST), WITHIN-LINEAGE#
	#   Dating is valid only inside one clonal lineage, so this #
	#   scopes to PA445's ST. Steps:                            #
	#     8a select focal-ST taxa + tip dates                   #
	#     8b recombination-free SNP alignment (snippy + Gubbins)#
	#     8c temporal-signal check (TreeTime) - GATE            #
	#     8d generate BEAST XML                                 #
	#     8e run BEAST                                          #
	#     8f TreeAnnotator -> MCC time tree                     #
	#############################################################

	echo -e "\e[33m ========================================== \e[0m"
	echo -e "\e[33m PHASE 8: TIME-CALIBRATED PHYLOGENY (BEAST) \e[0m"
	echo -e "\e[33m ========================================== \e[0m"

	# --- 8a. focal-ST taxa + tip dates -------------------------------------
	echo -e "\e[31m ===================================== \e[0m"
	echo -e "\e[31m Selecting focal-ST taxa and tip dates \e[0m"
	echo -e "\e[31m ===================================== \e[0m"

	BEAST_TAXA="${BEAST_PA}/focal_taxa.txt"
	BEAST_DATES="${BEAST_PA}/tip_dates.tsv"
	
	conda run -n beast_env python3 "${prep_beast_inputs}" \
	--mlst "${MLST_PA}/all.mlst.tsv" \
	--metadata "${REFERENCES_PA}/download_logs/accessions.tsv" \
	--focal "PA445" \
	--focal-date "${PA445_date}" \
	${focal_st:+--focal-st "${focal_st}"} \
	--min-strains 5 \
	--out-strains "${BEAST_TAXA}" \
	--out-dates "${BEAST_DATES}" \
	--out-report "${BEAST_PA}/tip_date_report.tsv"

	# --- 8a-guard: PA445 must carry a tip date, or it is silently dropped from the
	# alignment-vs-dates intersection in make_beast_xml (empty --focal-date is the
	# usual cause). Fail loudly here rather than after a multi-hour BEAST run.
	if ! awk -F'\t' 'NR>1{print $1}' "${BEAST_DATES}" | grep -qx "PA445"; then
		echo -e "\e[32m   ERROR: PA445 absent from ${BEAST_DATES} - check \${PA445_date}/--focal-date and tip_date_report.tsv. \e[0m"
		exit 1
	fi
	if ! grep -qx "PA445" "${BEAST_TAXA}"; then
		echo -e "\e[32m   ERROR: PA445 absent from ${BEAST_TAXA}. \e[0m"
		exit 1
	fi
	echo -e "\e[32m   PA445 present in tip dates and taxa set. \e[0m"

	# --- 8b. recombination-free SNP alignment ------------------------------
	# Map each focal isolate to PA445 (complete reference) and build a core
	# alignment, then remove recombination with Gubbins.
	echo -e "\e[31m ===================================== \e[0m"
	echo -e "\e[31m Snippy mapping vs PA445, then Gubbins \e[0m"
	echo -e "\e[31m ===================================== \e[0m"
	
	REF_FASTA="${REFERENCES_PA}/PA445.fasta"
	
	while IFS= read -r strain_name; do
		[[ -z "${strain_name}" ]] && continue
		[[ "${strain_name}" == "PA445" ]] && continue
		
		qfasta="${REFERENCES_PA}/${strain_name}.fasta"
		[[ -f "${qfasta}" ]] || { echo -e "\e[32m WARNING: ${qfasta} missing - skipping \e[0m"; continue; }
		
		conda run -n snippy_env snippy \
		--outdir "${SNIPPY_PA}/${strain_name}" \
		--ref "${REF_FASTA}" \
		--ctgs "${qfasta}" \
		--cpus "${threads}" --force
	done < "${BEAST_TAXA}"
	conda run -n snippy_env --no-capture-output bash -c \
	"cd '${SNIPPY_PA}' && snippy-core --ref '${REF_FASTA}' --prefix core */"

	# snippy-core always labels the --ref genome's own row literally "Reference" in its output alignments, regardless of the reference FASTA's filename.
	# PA445 (skipped above as the mapping target, not a query) would therefore be silently invisible to every downstream step that looks it up by name
	# (Gubbins, TreeTime, and make_beast_xml's date lookup) - it does NOT get dropped from the data, just from anything matching on "PA445" as a string.
	# Rename that one record back to PA445 immediately so the correct name propagates through Gubbins, TreeTime, and the BEAST XML.
	CORE_ALN="${SNIPPY_PA}/core.full.aln"
	if ! grep -q '^>Reference$' "${CORE_ALN}"; then
		echo -e "\e[32m   ERROR: expected a '>Reference' record in ${CORE_ALN} (snippy-core's usual label for the --ref genome) but none was found - snippy-core's output format may differ in your version; inspect ${CORE_ALN} before trusting the dating run. \e[0m"
		exit 1
	fi
	awk '/^>Reference$/{print ">PA445"; next} {print}' "${CORE_ALN}" > "${CORE_ALN}.tmp" \
	&& mv "${CORE_ALN}.tmp" "${CORE_ALN}"
	echo -e "\e[32m   Renamed snippy-core's 'Reference' row -> 'PA445' in core.full.aln \e[0m"

	# Gubbins on the whole-genome core alignment (needs genome coordinates)
	conda run -n beast_env  --no-capture-output bash -c \
	"cd '${GUBBINS_PA}' && run_gubbins.py --prefix gubbins --threads ${threads} '${SNIPPY_PA}/core.full.aln'"

	# recombination-free SNP alignment + constant-site counts (for ascertainment)
	SNP_ALN="${BEAST_PA}/core.recombfree.snps.fasta"
	conda run -n beast_env snp-sites -c \
	"${GUBBINS_PA}/gubbins.filtered_polymorphic_sites.fasta" > "${SNP_ALN}"
	
	conda run -n beast_env snp-sites -C \
	"${SNIPPY_PA}/core.full.aln" > "${BEAST_PA}/constant_site_counts.txt" 2>/dev/null || true
	
	echo -e "\e[32m   SNP alignment: ${SNP_ALN}; constant counts: ${BEAST_PA}/constant_site_counts.txt \e[0m"

	# --- 8c. temporal-signal gate (root-to-tip via TreeTime) ---------------
	echo -e "\e[31m ===================================================== \e[0m"
	echo -e "\e[31m Temporal-signal check (inspect before trusting dates) \e[0m"
	echo -e "\e[31m ===================================================== \e[0m"
	
	# TreeTime wants a name,date CSV
	awk -F'\t' 'NR>1{print $1","$2}' "${BEAST_DATES}" > "${BEAST_PA}/treetime_dates.csv"
	sed -i '1i name,date' "${BEAST_PA}/treetime_dates.csv"
	
	conda run -n beast_env  --no-capture-output bash -c \
	"cd '${BEAST_PA}' && treetime clock --tree '${GUBBINS_PA}/gubbins.final_tree.tre' \
	 --dates treetime_dates.csv --aln '${SNP_ALN}' --outdir treetime_clock" \
	 || echo -e "\e[32m   TreeTime clock check did not complete - inspect inputs \e[0m"
	echo -e "\e[32m   >>> Review ${BEAST_PA}/treetime_clock/ : if root-to-tip R^2 is low or the clock rate is negative, the lineage is NOT datable - do not trust BEAST. \e[0m"

	# --- 8d. generate BEAST XML --------------------------------------------
	echo -e "\e[31m ==================================================== \e[0m"
	echo -e "\e[31m Generating BEAST XML (${beast_clock} clock, ${beast_tree_prior} prior) \e[0m"
	echo -e "\e[31m ==================================================== \e[0m"

	BEAST_XML="${BEAST_PA}/PA445.lineage.beast.xml"
	
	conda run -n beast_env python3 "${make_beast_xml}" \
	--alignment "${SNP_ALN}" \
	--dates "${BEAST_DATES}" \
	--out "${BEAST_XML}" \
	--clock "${beast_clock}" \
	--tree-prior "${beast_tree_prior}" \
	--chain-length "${beast_chain_length}" \
	--log-every "${beast_log_every}" \
	--prefix "PA445.lineage"

	# Confirm PA445 is an actual taxon in the XML (taxon="PA445"), not just present
	# as the file prefix (PA445.lineage.log/.trees), which a bare grep would match.
	if ! grep -q 'taxon="PA445"' "${BEAST_XML}"; then
		echo -e "\e[31m   ERROR: PA445 is not a taxon in ${BEAST_XML} - inspect the make_beast_xml stderr 'dropped' warning. \e[0m"
		exit 1
	fi
	echo -e "\e[32m   PA445 confirmed as a taxon in ${BEAST_XML}. \e[0m"

	# --- 8e. run BEAST ------------------------------------------------------
	if [ "${run_beast}" = "true" ]; then
		echo -e "\e[31m =============================================================== \e[0m"
		echo -e "\e[31m Running BEAST (this is long; logs/trees land in ${BEAST_PA}) \e[0m"
		echo -e "\e[31m =============================================================== \e[0m"
			
		conda run -n beast_env --no-capture-output bash -c \
		"cd '${BEAST_PA}' && beast "-threads" ${threads} "-overwrite" '${BEAST_XML}'"

		# --- 8f. MCC time tree ---------------------------------------------
		echo -e "\e[32m TreeAnnotator -> MCC tree \e[0m"
		conda run -n beast_env --no-capture-output bash -c \
		"cd '${BEAST_PA}' && treeannotator -burnin ${beast_burnin_pct} -height median \
		 PA445.lineage.trees PA445.lineage.mcc.tree"
		echo -e "\e[32m Time-calibrated MCC tree: ${BEAST_PA}/PA445.lineage.mcc.tree \e[0m"
		echo -e "\e[32m >>> BEFORE USING: open PA445.lineage.log in Tracer and confirm ESS>200 for \e[0m"
		echo -e "\e[32m >>> all parameters. If not converged, raise beast_chain_length and re-run. \e[0m"
	else
		echo -e "\e[32m 8e: run_beast=false - XML ready at ${BEAST_XML} (verify in BEAUti, then run BEAST). \e[0m"
	fi

	echo ""
echo -e "\e[32m ======================================================================== \e[0m"
echo -e "\e[32m COMPARATIVE GENOMICS COMPLETE - $(date) \e[0m"
echo -e "\e[32m ======================================================================== \e[0m"
echo ""
echo -e "\e[32m -- FULL DATASET ------------------------------------------------------ \e[0m"
echo -e "\e[32m  Genomes downloaded : ${REFERENCES_PA}/ \e[0m"
echo -e "\e[32m  MLST results       : ${MLST_PA}/all.mlst.tsv \e[0m"
echo -e "\e[32m  Pasty serotypes    : ${PASTY_PA}/PA_serotypes.tsv \e[0m"
echo -e "\e[32m  RGI (all)          : ${RGI_PA}/ \e[0m"
echo -e "\e[32m  Prokka annotations : ${PROKKA_PA}/ \e[0m"
echo -e "\e[32m  Pan-genome         : ${PANAROO_PA}/results/ \e[0m"
echo ""
echo -e "\e[32m -- NDM-1 CONFIRMATION ----------------------------------------------- \e[0m"
echo -e "\e[32m  True NDM-1 strains : ${NDM1_STRAINS} \e[0m"
echo -e "\e[32m  Reconcile output   : ${CLINKER_PA}/truepos_gbk/ \e[0m"
echo ""
echo -e "\e[32m -- SUBSAMPLED DATASET ----------------------------------------------- \e[0m"
echo -e "\e[32m  Strain list        : ${SUBSAMPLE_LIST} \e[0m"
echo -e "\e[32m  FASTAs             : ${SUBSAMPLE_FASTA}/ \e[0m"
echo -e "\e[32m  RGI heatmap        : ${RGI_SUB}/PA445.AMR_heatmap_subsampled \e[0m"
echo -e "\e[32m  Phylogeny tree     : ${IQTREE_SUB}/PA445.phylogeny_subsampled.treefile \e[0m"
echo ""
echo -e "\e[32m -- CHROMOSOMAL-NDM-1 CONTROLLED VIEW (PHASE 7) ---------------------- \e[0m"
echo -e "\e[32m  blaNDM-1 location  : ${SUBSUB_PA}/ndm1_location_report.tsv \e[0m"
echo -e "\e[32m  Chromosomal NDM-1  : ${SUBSUB_PA}/chromosomal_ndm1_strains.txt \e[0m"
echo -e "\e[32m  Controlled set     : ${SUBSUB_PA}/subsubsampled_strains.txt \e[0m"
echo -e "\e[32m  RGI heatmap        : ${RGI_SUBSUB}/PA445.AMR_heatmap_chromosomal_ndm1 \e[0m"
echo -e "\e[32m  Phylogeny tree     : ${IQTREE_SUBSUB}/PA445.phylogeny_chromosomal_ndm1.treefile \e[0m"
echo ""
echo -e "\e[32m -- TIME-CALIBRATED PHYLOGENY (PHASE 8, within-lineage) ------------- \e[0m"
echo -e "\e[32m  Tip dates / report : ${BEAST_PA}/tip_dates.tsv , tip_date_report.tsv \e[0m"
echo -e "\e[32m  Recomb-free SNPs   : ${BEAST_PA}/core.recombfree.snps.fasta \e[0m"
echo -e "\e[32m  Temporal signal    : ${BEAST_PA}/treetime_clock/ \e[0m"
echo -e "\e[32m  BEAST XML          : ${BEAST_PA}/PA445.lineage.beast.xml \e[0m"
echo -e "\e[32m  MCC time tree      : ${BEAST_PA}/PA445.lineage.mcc.tree (after convergence check) \e[0m"
echo ""
echo -e "\e[32m -- NDM-1 STRAINS ONLY ---------------------------------------------- \e[0m"
echo -e "\e[32m  GBK region files   : ${CLINKER_PA}/gbk_inputs/ \e[0m"
echo -e "\e[32m  Synteny plot       : ${CLINKER_PA}/blaNDM1.clinker.html \e[0m"
echo -e "\e[32m  Platon             : ${PLATON_PA}/ \e[0m"
echo -e "\e[32m  PlasmidFinder      : ${PLASMIDFINDER_PA}/ \e[0m"
echo -e "\e[32m  MOB-suite          : ${MOBSUITE_PA}/ \e[0m"
echo ""
echo -e "\e[32m -- SUGGESTED NEXT STEPS -------------------------------------------- \e[0m"
echo -e "\e[32m  1. MLST: identify PA445 sequence type vs global NDM-1 strains \e[0m"
echo -e "\e[32m  2. RGI heatmap: compare AMR profiles in subsampled set \e[0m"
echo -e "\e[32m  3. Phylogeny: load treefile in FigTree/iTOL; annotate NDM-1 carriers \e[0m"
echo -e "\e[32m  4. Synteny: open blaNDM1.clinker.html in browser \e[0m"
echo -e "\e[32m  5. Plasmid: check MOB-suite chromosome.fasta for chromosomal blaNDM-1 \e[0m"
echo -e "\e[32m     (consistent with PA445 chromosomal integration) \e[0m"
echo -e "\e[32m  Full log: ${LOG} \e[0m"
echo -e "\e[32m ======================================================================== \e[0m"