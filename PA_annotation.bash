#!/bin/bash
# Script for annotating P. aeruginosa PA445 hybrid assembly and investigating blaNDM-1 duplication mechanism
# Building No.3
set -euo pipefail
################
#GLOBAL SETTING#
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/Paeruginosa"
ANNOTATION_PA445="${SAMPLE_PATH}/annotation_Paeruginosa"

#paths to databases
blast_db="${REF_PATH}/miniconda3/envs/BPannotate/db/blast/mlst.fa"
pubmlst="${REF_PATH}/miniconda3/envs/BPannotate/db/pubmlst"
bakta_db="${REF_PATH}/bakta_db/db"

#Create workflow directories
CHECKM_PA445="${ANNOTATION_PA445}/checkm"
CHECKM_INPUTS="${ANNOTATION_PA445}/checkm/checkm_inputs"
BAKTA_PA445="${ANNOTATION_PA445}/bakta"
GENE_PA445="${ANNOTATION_PA445}/genes"
MGE_PA445="${ANNOTATION_PA445}/mge"
MAFFT_PA445="${ANNOTATION_PA445}/mafft"
BLAST_PA445="${ANNOTATION_PA445}/blast"
BLAST_PA445_INPUTS="${ANNOTATION_PA445}/blast/inputs"
EXT_PATH="${ANNOTATION_PA445}/mge/extended_analysis"
RECOMBINATION_PA445="${ANNOTATION_PA445}/recombination"
PHISPY_PA445="${ANNOTATION_PA445}/phispy"
SYNTENY_PA445="${ANNOTATION_PA445}/synteny"

mkdir -p "${SAMPLE_PATH}/annotation_Paeruginosa"
mkdir -p "${ANNOTATION_PA445}/checkm"
mkdir -p "${ANNOTATION_PA445}/checkm/checkm_inputs"
mkdir -p "${CHECKM_PA445}/temp"
mkdir -p "${ANNOTATION_PA445}/bakta"
mkdir -p "${ANNOTATION_PA445}/genes"
mkdir -p "${ANNOTATION_PA445}/mge"
mkdir -p "${ANNOTATION_PA445}/mafft"
mkdir -p "${ANNOTATION_PA445}/blast"
mkdir -p "${ANNOTATION_PA445}/blast/inputs"
mkdir -p "${EXT_PATH}/ice"
mkdir -p "${EXT_PATH}/is91_boundaries"
mkdir -p "${ANNOTATION_PA445}/recombination"
mkdir -p "${ANNOTATION_PA445}/phispy"
mkdir -p "${ANNOTATION_PA445}/synteny/gbk_inputs"
#mkdir -p "${EXT_PATH}/logs"
#mkdir -p "${EXT_PATH}/extended_region"

#Other settings
gene_name="blaNDM-1"
flanking_bp=20000
threads=16

#ICE_B81 accession - NDM-P ICE_B81 from P. aeruginosa
#Update this accession if a more specific match is found from NCBI BLAST
ICE_B81_ACCESSION="OQ806931"

#Assembly paths
filtered_assembly="${SAMPLE_PATH}/445.hybrid_assembly_filter500bp.fasta"
filtered_assembly_renamed="${SAMPLE_PATH}/445.hybrid_assembly_filter500bp.renamed.fasta"

#Generate logs
LOG="${ANNOTATION_PA445}/PA_annotation.pipeline.log"
exec > >(tee -a "${LOG}") 2>&1

echo "========================================================"
echo " PA445 Annotation Pipeline - Started: $(date)"
echo "========================================================"

	#############################################
	#RENAME FASTA HEADER FOR CONSISTENT contig_1#
	#############################################

	#BUG FIX 2: Rename done after directories created, before any tools run
	#All tools use filtered_assembly_renamed so all outputs share contig_1 ID
	sed 's/^>1 .*/\>contig_1/' \
	"${filtered_assembly}" \
	> "${filtered_assembly_renamed}"

	echo -e "\e[32m Verification of renamed FASTA header: \e[0m"
	grep ">" "${filtered_assembly_renamed}"

	###########################
	#CHECK GENOME COMPLETENESS#
	###########################

	echo -e "\e[31m ============= \e[0m"
	echo -e "\e[31m CHECKM: PA445 \e[0m"
	echo -e "\e[31m ============= \e[0m"

	#BUG FIX 3: CheckM uses renamed assembly for contig_1 consistency
	cp "${filtered_assembly_renamed}" "${CHECKM_INPUTS}/"

	#Checking genome completeness
	conda run -n BPstructure checkm lineage_wf \
	-t ${threads} \
	--reduced_tree \
	--pplacer_threads 1 \
	-x fasta \
	--tmpdir "${CHECKM_PA445}/temp" \
	--tab_table \
	"${CHECKM_INPUTS}" \
	"${CHECKM_PA445}"

	#CheckM quality assessment and summary
	if [ -f "${CHECKM_PA445}/lineage.ms" ]; then
		conda run -n BPstructure checkm qa \
		"${CHECKM_PA445}/lineage.ms" \
		"${CHECKM_PA445}" \
		-o 2 \
		--tab_table \
		-f "${CHECKM_PA445}/PA445.quality.checkm.tsv"
		echo -e "\e[32m CheckM complete: ${CHECKM_PA445}/PA445.quality.checkm.tsv \e[0m"
	else
		echo -e "\e[32m WARNING: lineage.ms not found - skipping checkm qa \e[0m"
		echo -e "\e[32m Genome quality assessed by assembly statistics instead \e[0m"
	fi

	#######################
	#FUNCTIONAL ANNOTATION#
	#######################

	echo -e "\e[31m ============ \e[0m"
	echo -e "\e[31m BAKTA: PA445 \e[0m"
	echo -e "\e[31m ============ \e[0m"

	#BUG FIX 4: Bakta uses renamed assembly so GFF3/GBK outputs have contig_1
	conda run -n BPannotation bakta \
	--force \
	--threads ${threads} \
	--db "${bakta_db}" \
	--genus Pseudomonas \
	--species aeruginosa \
	--prefix "PA445.hybrid.filtered" \
	--output "${BAKTA_PA445}" \
	"${filtered_assembly_renamed}"

	######################
	#EXTRACT GENE REGIONS#
	######################

	#Determine blaNDM-1 gene copy locations from Bakta TSV
	echo -e "\e[32m Extracting ${gene_name} coordinates from Bakta annotation \e[0m"

	grep -i "${gene_name}" "${BAKTA_PA445}/PA445.hybrid.filtered.tsv" \
	> "${GENE_PA445}/${gene_name}.coordinates.txt" || true

	nhits=$(wc -l < "${GENE_PA445}/${gene_name}.coordinates.txt")

	#Log annotation hits
	{
	echo ""
	echo "Found ${nhits} annotation hit(s) for ${gene_name}:"
	cat "${GENE_PA445}/${gene_name}.coordinates.txt"
	} >> "${GENE_PA445}/${gene_name}.coordinates.log" 2>&1

	echo -e "\e[32m Found ${nhits} hit(s) for ${gene_name} - see ${GENE_PA445}/${gene_name}.coordinates.log \e[0m"

	if [ "${nhits}" -eq 0 ]; then
		echo -e "\e[32m WARNING: No hits found for ${gene_name} in Bakta TSV \e[0m"
		echo -e "\e[32m Check the gene name or inspect ${BAKTA_PA445}/PA445.hybrid.filtered.tsv manually \e[0m"
		exit 1
	fi

	echo -e "\e[31m ============= \e[0m"
	echo -e "\e[31m SEQKIT: PA445 \e[0m"
	echo -e "\e[31m ============= \e[0m"

	#Step 1: Extract each gene copy sequence individually
	copy=1

while IFS=$'\t' read -r contig type start end strand locus gene product rest;
do
	{
	echo ""
	echo "--- Copy ${copy} ---"
	echo "Contig  : ${contig}"
	echo "Start   : ${start}"
	echo "End     : ${end}"
	echo "Strand  : ${strand}"
	echo "Locus   : ${locus}"
	echo "Product : ${product}"
	} >> "${GENE_PA445}/${gene_name}.coordinates.log" 2>&1

	echo -e "\e[32m Extracting copy ${copy}: ${contig}:${start}-${end} (${strand}) \e[0m"

	conda run -n assembly seqkit subseq \
	-r "${start}:${end}" \
	"${filtered_assembly_renamed}" \
	> "${GENE_PA445}/${gene_name}_copy${copy}.fasta"

	if grep -qi "unknown flag\|error\|Error" \
	"${GENE_PA445}/${gene_name}_copy${copy}.fasta" 2>/dev/null; then
		echo -e "\e[32m ERROR: seqkit failed for copy ${copy} - contents of output file: \e[0m"
		cat "${GENE_PA445}/${gene_name}_copy${copy}.fasta"
		exit 1
	fi

	echo -e "\e[32m Saved: ${GENE_PA445}/${gene_name}_copy${copy}.fasta \e[0m"

	copy=$(( copy + 1 ))

done < "${GENE_PA445}/${gene_name}.coordinates.txt"

	#Step 2: Extract combined region spanning both copies + flanks
	echo ""
	echo -e "\e[32m Extracting combined region spanning both ${gene_name} copies... \e[0m"

	region_start=$(awk -F'\t' 'NR==1 {print $3}' "${GENE_PA445}/${gene_name}.coordinates.txt")
	region_end=$(awk -F'\t' 'END {print $4}' "${GENE_PA445}/${gene_name}.coordinates.txt")
	region_contig=$(awk -F'\t' 'NR==1 {print $1}' "${GENE_PA445}/${gene_name}.coordinates.txt")

	combined_start=$(( region_start - flanking_bp > 0 ? region_start - flanking_bp : 1 ))
	combined_end=$(( region_end + flanking_bp ))

	echo -e "\e[32m Contig         : ${region_contig} \e[0m"
	echo -e "\e[32m Copy 1 start   : ${region_start} \e[0m"
	echo -e "\e[32m Copy 2 end     : ${region_end} \e[0m"
	echo -e "\e[32m Combined region: ${combined_start}-${combined_end} (with +/-${flanking_bp} bp flanks) \e[0m"

	#Log combined region coordinates
	{
	echo ""
	echo "Combined region: ${region_contig}:${combined_start}-${combined_end}"
	echo "Spans both copies + ${flanking_bp} bp flanks on each side"
	} >> "${GENE_PA445}/${gene_name}.coordinates.log" 2>&1

	conda run -n assembly seqkit subseq \
	-r "${combined_start}:${combined_end}" \
	"${filtered_assembly_renamed}" \
	> "${GENE_PA445}/${gene_name}.combined_region.fasta"

	if grep -qi "unknown flag\|error\|Error" \
	"${GENE_PA445}/${gene_name}.combined_region.fasta" 2>/dev/null; then
		echo -e "\e[32m ERROR: seqkit failed for combined region - contents of output file: \e[0m"
		cat "${GENE_PA445}/${gene_name}.combined_region.fasta"
		exit 1
	fi

	echo -e "\e[32m Saved: ${GENE_PA445}/${gene_name}.combined_region.fasta \e[0m"
	echo ""
	echo -e "\e[32m Combined region structure: \e[0m"
	echo -e "\e[32m |--[flank ${flanking_bp}bp]--[copy1]--[intergenic]--[copy2]--[flank ${flanking_bp}bp]--| \e[0m"

	##############################
	#IDENTIFY INSERTION SEQUENCES#
	##############################

	flank_combined="${GENE_PA445}/${gene_name}.combined_region.fasta"

	#--- ISEScan: combined region (focused analysis) ---
	echo -e "\e[31m ================================ \e[0m"
	echo -e "\e[31m ISESCAN: PA445 (combined region) \e[0m"
	echo -e "\e[31m ================================ \e[0m"

	conda run -n recombination isescan.py \
	--nthread ${threads} \
	--seqfile "${flank_combined}" \
	--output "${MGE_PA445}/${gene_name}.combined_region.isescan"

	#--- IntegronFinder: combined region (focused analysis) ---
	echo -e "\e[31m ======================================== \e[0m"
	echo -e "\e[31m INTEGRON FINDER: PA445 (combined region) \e[0m"
	echo -e "\e[31m ======================================== \e[0m"

	conda run -n recombination integron_finder \
	--gbk --pdf \
	--linear \
	--outdir "${MGE_PA445}/${gene_name}.combined_region.integron" \
	"${flank_combined}"

	#Rename IntegronFinder combined region outputs
	INTEGRON_RESULTS="${MGE_PA445}/${gene_name}.combined_region.integron/Results_Integron_Finder_${gene_name}.combined_region"

	echo -e "\e[32m Renaming IntegronFinder combined region outputs... \e[0m"

	#BUG FIX 5: check and move contig_1.gbk (not 1.gbk) since renamed FASTA used
	if [ -f "${INTEGRON_RESULTS}/contig_1.gbk" ]; then
		mv "${INTEGRON_RESULTS}/contig_1.gbk" \
		"${INTEGRON_RESULTS}/${gene_name}.combined_region.integron.gbk"
		echo -e "\e[32m Renamed: contig_1.gbk to ${gene_name}.combined_region.integron.gbk \e[0m"
	else
		echo -e "\e[32m WARNING: contig_1.gbk not found - skipping GBK rename \e[0m"
		echo -e "\e[32m Expected if no integrons found in combined region \e[0m"
	fi

	pdf_count=0

for pdf_file in "${INTEGRON_RESULTS}"/contig_1_*.pdf;
do
	[ -f "${pdf_file}" ] || continue
	#BUG FIX 6: use rev to correctly extract integron number from contig_1_N.pdf
	integron_num=$(basename "${pdf_file}" .pdf | rev | cut -d'_' -f1 | rev)
	new_name="${gene_name}.combined_region.integron_${integron_num}.pdf"
	mv "${pdf_file}" "${INTEGRON_RESULTS}/${new_name}"
	echo -e "\e[32m Renamed: $(basename ${pdf_file}) to ${new_name} \e[0m"
	pdf_count=$(( pdf_count + 1 ))
done

	echo -e "\e[32m Total PDFs renamed: ${pdf_count} \e[0m"
	echo -e "\e[32m IntegronFinder combined region output files: \e[0m"
	ls -lh "${INTEGRON_RESULTS}/"

	#--- ISEScan: full genome (supportive analysis) ---
	echo -e "\e[31m ============================ \e[0m"
	echo -e "\e[31m ISESCAN: PA445 (full genome) \e[0m"
	echo -e "\e[31m ============================ \e[0m"

	#BUG FIX 7: ISEScan full genome uses renamed assembly for contig_1 consistency
	conda run -n recombination isescan.py \
	--nthread ${threads} \
	--seqfile "${filtered_assembly_renamed}" \
	--output "${MGE_PA445}/PA445.full_genome.isescan"

	#--- IntegronFinder: full genome (supportive analysis) ---
	echo -e "\e[31m ==================================== \e[0m"
	echo -e "\e[31m INTEGRON FINDER: PA445 (full genome) \e[0m"
	echo -e "\e[31m ==================================== \e[0m"

	conda run -n recombination integron_finder \
	--gbk --pdf \
	--circ \
	--outdir "${MGE_PA445}/PA445.full_genome.integron" \
	"${filtered_assembly_renamed}"

	#Rename IntegronFinder full genome outputs
	#BUG FIX 8: Results dir named after renamed FASTA basename
	INTEGRON_RESULTS="${MGE_PA445}/PA445.full_genome.integron/Results_Integron_Finder_445.hybrid_assembly_filter500bp.renamed"

	echo -e "\e[32m Renaming IntegronFinder full genome outputs... \e[0m"

	#BUG FIX 9: check and move contig_1.gbk since renamed FASTA used
	if [ -f "${INTEGRON_RESULTS}/contig_1.gbk" ]; then
		mv "${INTEGRON_RESULTS}/contig_1.gbk" \
		"${INTEGRON_RESULTS}/PA445.hybrid.filtered.integrons.gbk"
		echo -e "\e[32m Renamed: contig_1.gbk to PA445.hybrid.filtered.integrons.gbk \e[0m"
	else
		echo -e "\e[32m WARNING: contig_1.gbk not found - skipping GBK rename \e[0m"
	fi

	pdf_count=0

for pdf_file in "${INTEGRON_RESULTS}"/contig_1_*.pdf;
do
	[ -f "${pdf_file}" ] || continue
	#BUG FIX 10: use rev to correctly extract integron number from contig_1_N.pdf
	integron_num=$(basename "${pdf_file}" .pdf | rev | cut -d'_' -f1 | rev)
	new_name="PA445.hybrid.filtered.integron_${integron_num}.pdf"
	mv "${pdf_file}" "${INTEGRON_RESULTS}/${new_name}"
	echo -e "\e[32m Renamed: $(basename ${pdf_file}) to ${new_name} \e[0m"
	pdf_count=$(( pdf_count + 1 ))
done

	echo -e "\e[32m Total PDFs renamed: ${pdf_count} \e[0m"
	echo -e "\e[32m IntegronFinder full genome output files: \e[0m"
	ls -lh "${INTEGRON_RESULTS}/"

	##############################
	#MULTIPLE SEQUENCES ALIGNMENT#
	##############################

	echo -e "\e[31m ============================== \e[0m"
	echo -e "\e[31m MAFFT: PA445 (blaNDM-1 COPIES) \e[0m"
	echo -e "\e[31m ============================== \e[0m"

	#Combine all copies into one file as MAFFT input
	cat "${GENE_PA445}/${gene_name}_copy"*.fasta \
	> "${MAFFT_PA445}/${gene_name}.all_copies.fasta"

	#Multiple sequence alignment (fasta output)
	conda run -n phylogeny mafft --auto \
	"${MAFFT_PA445}/${gene_name}.all_copies.fasta" \
	> "${MAFFT_PA445}/${gene_name}.aligned.fasta"

	echo -e "\e[32m Alignment saved: ${MAFFT_PA445}/${gene_name}.aligned.fasta \e[0m"

	#Multiple sequence alignment (ClustalW output)
	conda run -n phylogeny mafft --auto \
	--clustalout \
	"${MAFFT_PA445}/${gene_name}.all_copies.fasta" \
	> "${MAFFT_PA445}/${gene_name}.aligned.aln"

	echo -e "\e[32m Alignment saved: ${MAFFT_PA445}/${gene_name}.aligned.aln \e[0m"

	conda run -n assembly seqkit fx2tab \
	--only-id --name --seq-hash \
	"${MAFFT_PA445}/${gene_name}.aligned.fasta" \
	> "${MAFFT_PA445}/${gene_name}.seq_hashes.txt"

echo "================================================"
echo " PA445 Extended Analysis - Started: $(date)"
echo "================================================"
echo " blaNDM-1 duplication mechanism investigation"
echo " Investigate IS (IS91 and new_269), ICE"
echo " Finding recombination and restriction sites"
echo " Finding prophage sequences"
echo "================================================"

	##################
	#NEW_269 ANALYSES#
	##################

	flank_combined="${GENE_PA445}/${gene_name}.combined_region.fasta"

	echo -e "\e[31m =================================== \e[0m"
	echo -e "\e[31m BLASTN: PA445 (new_269 TRANSPOSASE) \e[0m"
	echo -e "\e[31m =================================== \e[0m"

	#Extract new_269 sequence from combined region
	#new_269 left copy: using more robust approach by extracting based on genome coordinate
	conda run -n assembly seqkit subseq \
	-r 2056385:2059007 \
	"${filtered_assembly_renamed}" \
	> "${BLAST_PA445_INPUTS}/new_269_transposase.fasta"

	echo -e "\e[32m new_269 sequence extracted: ${BLAST_PA445_INPUTS}/new_269_transposase.fasta \e[0m"

	#Build local BLAST database from renamed assembly for contig_1 consistency
	conda run -n ncbi makeblastdb \
	-in "${filtered_assembly_renamed}" \
	-dbtype nucl \
	-out "${BLAST_PA445}/PA445_db"

	#BLAST new_269 against genome to find all copies
	conda run -n ncbi blastn \
	-query "${BLAST_PA445_INPUTS}/new_269_transposase.fasta" \
	-db "${BLAST_PA445}/PA445_db" \
	-outfmt "6 qseqid sseqid pident length sstart send evalue bitscore" \
	-perc_identity 80 \
	-out "${BLAST_PA445}/new_269_vs_genome.txt"

	echo -e "\e[32m BLAST results: ${BLAST_PA445}/new_269_vs_genome.txt \e[0m"
	echo ""
	cat "${BLAST_PA445}/new_269_vs_genome.txt"

	##################
	#ICE_B81 ANALYSES#
	##################

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m BLASTN: PA445 (ICE_B81) \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	echo ""
	echo -e "\e[32m BLASTN revealed flanking regions containing blaNDM-1 hits with ICE_B81 mobile element \e[0m"
	echo ""	

	#Step 1: Download ICE_B81 reference sequence from NCBI
	#Only run if accession has been updated from placeholder
	if [ "${ICE_B81_ACCESSION}" == "UPDATE_WITH_NCBI_ACCESSION" ]; then
		echo -e "\e[32m WARNING: ICE_B81_ACCESSION not set \e[0m"
		echo -e "\e[32m Please search NCBI for 'NDM-P ICE_B81 Pseudomonas aeruginosa' \e[0m"
		echo -e "\e[32m Update ICE_B81_ACCESSION at the top of this script and rerun \e[0m"
		echo -e "\e[32m Skipping ICE_B81 download and BLAST steps \e[0m"
	else
		echo -e "\e[32m Downloading ICE_B81 reference (${ICE_B81_ACCESSION})... \e[0m"

		conda run -n ncbi efetch \
		-db nuccore \
		-id "${ICE_B81_ACCESSION}" \
		-format fasta \
		> "${EXT_PATH}/ice/ICE_B81_reference.fasta"

		echo -e "\e[32m ICE_B81 reference downloaded \e[0m"

		#Step 2: BLAST ICE_B81 against PA445 genome
		echo -e "\e[32m BLASTing ICE_B81 against PA445 genome... \e[0m"

		conda run -n ncbi blastn \
		-query "${EXT_PATH}/ice/ICE_B81_reference.fasta" \
		-db "${BLAST_PA445}/PA445_db" \
		-outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
		-perc_identity 80 \
		-out "${EXT_PATH}/ice/ICE_B81_vs_PA445.txt"

		echo -e "\e[32m ICE_B81 BLAST results: \e[0m"
		echo "qseqid | sseqid | %id | length | qstart | qend | sstart | send | evalue"
		cat "${EXT_PATH}/ice/ICE_B81_vs_PA445.txt"

		#Step 3: Extract the matching region from PA445 for clinker comparison
		if [ -s "${EXT_PATH}/ice/ICE_B81_vs_PA445.txt" ]; then
			#Get the genomic coordinates of the ICE match
			ice_start=$(awk 'NR==1 {print $7}' "${EXT_PATH}/ice/ICE_B81_vs_PA445.txt")
			ice_end=$(awk 'NR==1 {print $8}' "${EXT_PATH}/ice/ICE_B81_vs_PA445.txt")

			echo -e "\e[32m ICE_B81 match in PA445: contig_1:${ice_start}-${ice_end} \e[0m"

			#Log ICE boundaries
			{
			echo "ICE_B81 match in PA445:"
			echo "  Start: ${ice_start}"
			echo "  End  : ${ice_end}"
			echo "  Span : $(( ice_end - ice_start )) bp"
			} >> "${EXT_PATH}/ice/ICE_boundaries.log"

			#Extract ICE region from PA445 for annotation
			conda run -n assembly seqkit subseq \
			-r "${ice_start}:${ice_end}" \
			"${filtered_assembly_renamed}" \
			> "${EXT_PATH}/ice/PA445_ICE_region.fasta"

			echo -e "\e[32m PA445 ICE region extracted: ${EXT_PATH}/ice/PA445_ICE_region.fasta \e[0m"
		else
			echo -e "\e[32m WARNING: No ICE_B81 BLAST hits found in PA445 \e[0m"
		fi
	fi

	#Step 4: Prepare ICEfinder submission file
	#ICEfinder web tool: https://bioinfo-mml.sjtu.edu.cn/ICEfinder/ICEfinder.html
	echo -e "\e[32m Preparing ICEfinder submission... \e[0m"

	cp "${filtered_assembly_renamed}" \
	"${EXT_PATH}/ice/PA445_for_ICEfinder_submission.fasta"

	echo ""
	echo -e "\e[32m ICEfinder submission instructions: \e[0m"
	echo -e "\e[32m   1. Go to https://bioinfo-mml.sjtu.edu.cn/ICEfinder/ICEfinder.html \e[0m"
	echo -e "\e[32m   2. Upload: ${EXT_PATH}/ice/PA445_for_ICEfinder_submission.fasta \e[0m"
	echo -e "\e[32m   3. Select organism type: Gram-negative \e[0m"
	echo -e "\e[32m   4. Download results and save to: ${EXT_PATH}/ice/icefinder_results/ \e[0m"

	###############
	#IS91 ANALYSES#
	###############

	echo -e "\e[31m ===================================== \e[0m"
	echo -e "\e[31m SEQKIT: IS91 (oriIS/terIS BOUNDARIES) \e[0m"
	echo -e "\e[31m ===================================== \e[0m"

	#Extract IS91 sequences from extended region using ISEScan coordinates
	#IS91(3) and IS91(4) positions from previous ISEScan analysis
	#Convert to genome coordinates and extract

	echo -e "\e[32m Extracting IS91 element sequences for boundary analysis \e[0m"

	#IS91(3) - full copy
	#Genome coords from previous analysis: 2059299-2061324
	conda run -n assembly seqkit subseq \
	-r 2059299:2061324 \
	"${filtered_assembly_renamed}" \
	> "${EXT_PATH}/is91_boundaries/IS91_3_full.fasta"

	#IS91(4) - truncated copy
	#Genome coords from previous analysis: 2066159-2068187
	conda run -n assembly seqkit subseq \
	-r 2066159:2068187 \
	"${filtered_assembly_renamed}" \
	> "${EXT_PATH}/is91_boundaries/IS91_4_truncated.fasta"

	echo -e "\e[31m =========== \e[0m"
	echo -e "\e[31m MAFFT: IS91 \e[0m"
	echo -e "\e[31m =========== \e[0m"

	echo -e "\e[32m  Align IS91(3) and IS91(4) to confirm 5'-truncation \e[0m"
	
	cat "${EXT_PATH}/is91_boundaries/IS91_3_full.fasta" \
	"${EXT_PATH}/is91_boundaries/IS91_4_truncated.fasta" \
	> "${EXT_PATH}/is91_boundaries/IS91_both_copies.fasta"

	#IS91 sequences alignment (fasta format)
	conda run -n phylogeny mafft --auto \
	"${EXT_PATH}/is91_boundaries/IS91_both_copies.fasta" \
	> "${EXT_PATH}/is91_boundaries/IS91_alignment.fasta"

	echo -e "\e[32m IS91 alignment saved: ${EXT_PATH}/is91_boundaries/IS91_alignment.fasta \e[0m"

	#IS91 sequences alignment (ClustalW format)
	conda run -n phylogeny mafft --auto \
	--clustalout \
	"${EXT_PATH}/is91_boundaries/IS91_both_copies.fasta" \
	> "${EXT_PATH}/is91_boundaries/IS91_alignment.aln"

	echo -e "\e[32m IS91 alignment saved: ${EXT_PATH}/is91_boundaries/IS91_alignment.aln \e[0m"
	echo -e "\e[32m Open in AliView to confirm 5-prime truncation of IS91(4) \e[0m"

	echo -e "\e[31m ============ \e[0m"
	echo -e "\e[31m BLASTN: IS91  \e[0m"
	echo -e "\e[31m ============ \e[0m"

	#BLAST IS91(3) against ISFinder database to confirm IS91 family
	#and identify oriIS/terIS boundary sequences
	conda run -n ncbi blastn \
	-query "${EXT_PATH}/is91_boundaries/IS91_3_full.fasta" \
	-db "${BLAST_PA445}/PA445_db" \
	-outfmt "6 qseqid sseqid pident length sstart send evalue bitscore" \
	-perc_identity 80 \
	-out "${EXT_PATH}/is91_boundaries/IS91_self_blast.txt"

	echo -e "\e[32m IS91 self-BLAST (copy number confirmation): \e[0m"
	cat "${EXT_PATH}/is91_boundaries/IS91_self_blast.txt"

	#############################
	#RECOMBINATION SITE ANALYSES#
	#############################

	echo -e "\e[31m ================================== \e[0m"
	echo -e "\e[31m EINVERTED: PA445 (combined region) \e[0m"
	echo -e "\e[31m ================================== \e[0m"

	#ext_region="${EXT_PATH}/extended_region/${gene_name}.extended_region.fasta"

	#--- Inverted repeats (einverted) ---
	echo -e "\e[32m Finding inverted repeats (einverted)... \e[0m"

	conda run -n recombination einverted \
	-sequence "${flank_combined}" \
	-gap 12 \
	-threshold 50 \
	-match 3 \
	-mismatch -4 \
	-outfile "${RECOMBINATION_PA445}/${gene_name}.combined_region.einverted.txt" \
	-outseq "${RECOMBINATION_PA445}/${gene_name}.combined_region.einverted.fasta"

	echo -e "\e[32m Inverted repeats: ${RECOMBINATION_PA445}/${gene_name}.combined_region.einverted.txt \e[0m"

	echo -e "\e[31m =================================== \e[0m"
	echo -e "\e[31m PALINDROME: PA445 (combined region) \e[0m"
	echo -e "\e[31m =================================== \e[0m"

	#--- Palindromes (palindrome) ---
	echo -e "\e[32m Finding palindromic restriction sites (palindrome)... \e[0m"

	conda run -n recombination palindrome \
	-auto \
	-sequence "${flank_combined}" \
	-minpallen 6 \
	-maxpallen 100 \
	-gaplimit 100 \
	-nummismatches 0 \
	-outfile "${RECOMBINATION_PA445}/${gene_name}.combined_region.palindrome.txt"

	echo -e "\e[32m Palindromes: ${RECOMBINATION_PA445}/${gene_name}.combined_region.palindrome.txt \e[0m"

	echo -e "\e[31m ================================ \e[0m"
	echo -e "\e[31m FUZZNUC: PA445 (combined region) \e[0m"
	echo -e "\e[31m ================================ \e[0m"

	#--- Direct repeats (dottup / fuzznuc) ---
	#Direct repeats flanking the duplicated cassette are evidence of
	#target site duplication from IS91 insertion
	echo -e "\e[32m Searching for target site duplications (fuzznuc)... \e[0m"

	#Extract the 20 bp flanking each IS91 boundary for target sites duplicatoin (TSD) search
	#IS91(3) left boundary: genome position 2059299
	#Extract 10 bp upstream and downstream of insertion site
	conda run -n recombination fuzznuc \
	-sequence "${flank_combined}" \
	-pattern "AAAAAAAAAA" \
	-complement Y \
	-outfile "${RECOMBINATION_PA445}/${gene_name}.fuzznuc.tsd_candidates.txt" || true

	echo -e "\e[32m TSD candidates: ${RECOMBINATION_PA445}/${gene_name}.fuzznuc.tsd_candidates.txt \e[0m"

	#--- att site search using IS91-specific terminal sequences ---
	#IS91 uses specific 9-bp AT-rich sequences at oriIS and terIS
	#Search for the consensus oriIS sequence in the extended region

	echo -e "\e[32m Searching for IS91 oriIS consensus sequence... \e[0m"

	conda run -n recombination fuzznuc \
	-sequence "${flank_combined}" \
	-pattern "ATAATTTTT" \
	-complement Y \
	-outfile "${RECOMBINATION_PA445}/PA445.combined_region.fuzznuc.IS91_oriIS_candidates.txt" || true

	echo -e "\e[32m oriIS candidates: ${RECOMBINATION_PA445}/PA445.combined_region.fuzznuc.IS91_oriIS_candidates.txt \e[0m"

	###################
	#PROPHAGE ANALYSES#
	###################

	echo -e "\e[31m ============= \e[0m"
	echo -e "\e[31m PHISPY: PA445 \e[0m"
	echo -e "\e[31m ============= \e[0m"

	conda run -n recombination PhiSpy.py \
	"${BAKTA_PA445}/PA445.hybrid.filtered.gbff" \
	-o "${PHISPY_PA445}/PA445.full.prophages" \
	--output_choice 7

	#########
	#SYNTENY#
	#########

	echo -e "\e[31m ====================== \e[0m"
	echo -e "\e[31m CLINKER: PA445 SYNTENY \e[0m"
	echo -e "\e[31m ====================== \e[0m"

	#Copy PA445 Bakta GBK to synteny inputs
	cp "${BAKTA_PA445}/PA445.hybrid.filtered.gbff" \
	"${SYNTENY_PA445}/gbk_inputs/PA445.gbk"

	#Extract just the blaNDM-1 extended region from PA445 GBK
	#for a focused synteny comparison matching supervisor's figure style
	#This requires extracting the relevant GBK record for just the region

conda run -n assembly python3 << PYEOF
from Bio import SeqIO

#Load full Bakta GBK
record = SeqIO.read("${BAKTA_PA445}/PA445.hybrid.filtered.gbk", "genbank")

#Extract extended region (0-based for BioPython)
region = record[${combined_start}-1:${combined_end}]
region.id = "PA445_blaNDM1_region"
region.name = "PA445_blaNDM1_region"
region.description = "PA445 blaNDM-1 extended region ${combined_start}-${combined_end}"

SeqIO.write(region, "${SYNTENY_PA445}/gbk_inputs/PA445_blaNDM1_region.gbk", "genbank")
print("Extracted region GBK saved")
print(f"Region size: {len(region)} bp")
print(f"Features: {len(region.features)}")
PYEOF

	echo -e "\e[32m PA445 blaNDM-1 region GBK: ${SYNTENY_PA445}/gbk_inputs/PA445_blaNDM1_region.gbk \e[0m"

	#Check how many GBK files are available for comparison
	gbk_count=$(ls "${SYNTENY_PA445}/gbk_inputs/"*.gbk 2>/dev/null | wc -l)
	echo -e "\e[32m GBK files available for synteny: ${gbk_count} \e[0m"
	ls -lh "${SYNTENY_PA445}/gbk_inputs/"

	echo ""
	echo -e "\e[32m To add reference strain comparison to clinker: \e[0m"
	echo -e "\e[32m   1. Identify blaNDM-1 carrying strains from AMR results \e[0m"
	echo -e "\e[32m   2. Download their genomes \e[0m"
	echo -e "\e[32m   3. Extract blaNDM-1 region using same seqkit approach \e[0m"
	echo -e "\e[32m   4. Annotate with Prokka to get GBK \e[0m"
	echo -e "\e[32m   5. Copy GBK to: ${SYNTENY_PA445}/gbk_inputs/ \e[0m"
	echo -e "\e[32m   6. Run clinker manually or rerun this script \e[0m"

	#Run clinker if more than 1 GBK present
	if [ "${gbk_count}" -gt 1 ]; then
		echo -e "\e[32m Running clinker on available GBK files... \e[0m"

		conda run -n phylogeny clinker \
		"${SYNTENY_PA445}/gbk_inputs/"*.gbk \
		--plot "${SYNTENY_PA445}/blaNDM1.combined_region.synteny.html" \
		--identity 0.3 \
		--genes

		echo -e "\e[32m Synteny plot: ${SYNTENY_PA445}/blaNDM1.combined_region.synteny.html \e[0m"
		echo -e "\e[32m Open in browser to view interactive figure \e[0m"
	else
		echo -e "\e[32m Only PA445 GBK available - add reference strain GBKs \e[0m"
		echo -e "\e[32m to ${SYNTENY_PA445}/gbk_inputs/ then rerun clinker \e[0m"
	fi

	################################
	#PREPARATION OF JBROWSE2 INPUTS#
	################################

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m JBROWSE2 INPUTS PREPARATION \e[0m"
	echo -e "\e[31m =========================== \e[0m"

	#Generate ISEScan coordinates in GFF format (combined region to genome coords)
	awk -F'\t' -v offset="${combined_start}" 'NR>1 {
	genome_start = $4 + offset - 1
	genome_end   = $5 + offset - 1
	counter++
	id = $2 "_" $3 "_" counter
	print "contig_1\tISEScan\tIS_element\t" genome_start "\t" genome_end \
	"\t.\t+\t.\tID=" id ";Name=" $2 "_" $3
	}' "${MGE_PA445}/${gene_name}.combined_region.isescan/genes/${gene_name}.combined_region.fasta.tsv" \
	> "${MGE_PA445}/${gene_name}.combined_region.isescan/${gene_name}.combined_region.isescan.genome_coords.unique.gff"

	echo -e "\e[32m ISEScan GFF generated: \e[0m"
	cat "${MGE_PA445}/${gene_name}.combined_region.isescan/${gene_name}.combined_region.isescan.genome_coords.unique.gff"

	#Convert IntegronFinder GBK to GFF3 (combined region)
	#genbank_to is installed in BPannotate env to avoid matplotlib conflict with integron_env
	INTEGRON_RESULTS="${MGE_PA445}/${gene_name}.combined_region.integron/Results_Integron_Finder_${gene_name}.combined_region"

	if [ -f "${INTEGRON_RESULTS}/${gene_name}.combined_region.integron.gbk" ]; then
		conda run -n ncbi genbank_to \
		-g "${INTEGRON_RESULTS}/${gene_name}.combined_region.integron.gbk" \
		--gff3 "${MGE_PA445}/${gene_name}.combined_region.integron/${gene_name}.combined_region.integron.gff3"
		echo -e "\e[32m GFF3 saved: ${MGE_PA445}/${gene_name}.combined_region.integron/${gene_name}.combined_region.integron.gff3 \e[0m"
	else
		echo -e "\e[32m WARNING: No GBK found for combined region - skipping genbank_to \e[0m"
		echo -e "\e[32m Expected if IntegronFinder found no integrons in combined region \e[0m"
	fi

	#Convert IntegronFinder GBK to GFF3 (full genome)
	INTEGRON_RESULTS="${MGE_PA445}/PA445.full_genome.integron/Results_Integron_Finder_445.hybrid_assembly_filter500bp.renamed"

	if [ -f "${INTEGRON_RESULTS}/PA445.hybrid.filtered.integrons.gbk" ]; then
		conda run -n ncbi genbank_to \
		-g "${INTEGRON_RESULTS}/PA445.hybrid.filtered.integrons.gbk" \
		--gff3 "${MGE_PA445}/PA445.full_genome.integron/PA445.hybrid.filtered.integrons.gff3"
		echo -e "\e[32m GFF3 saved: ${MGE_PA445}/PA445.full_genome.integron/PA445.hybrid.filtered.integrons.gff3 \e[0m"
	else
		echo -e "\e[32m WARNING: No GBK found for full genome - skipping genbank_to \e[0m"
	fi

	#ISEScan GFF already generated above with updated combined_start
	#echo -e "\e[32m ISEScan GFF (extended region, genome coords): \e[0m"
	#cat "${EXT_PATH}/extended_region/${gene_name}.extended_region.isescan.genome_coords.gff"

	#Generate recombination sites GFF for JBrowse2
	#Convert einverted output to GFF3
	echo -e "\e[32m Converting einverted output to GFF3... \e[0m"

	awk -v offset="${combined_start}" '
	/Score/ {
		#Parse score line: "contig_1: Score 54: 17/17 (100%) matches..."
		match($0, /Score ([0-9]+)/, arr)
		score = arr[1]
	}
	/^   [0-9]/ && score != "" {
		#Parse coordinate line
		match($0, /([0-9]+) [actgACTG]+ ([0-9]+)/, coords)
		if (coords[1] != "" && coords[2] != "") {
			start = coords[1] + offset - 1
			end   = coords[2] + offset - 1
			counter++
			print "contig_1\teinverted\tinverted_repeat\t" start "\t" end \
			"\t" score "\t.\t.\tID=IR_" counter ";Note=Score_" score
			score = ""
		}
	}' "${RECOMBINATION_PA445}/${gene_name}.combined_region.einverted.txt" \
	> "${RECOMBINATION_PA445}/${gene_name}.combined_region.einverted.genome_coords.gff"

	echo -e "\e[32m Inverted repeat GFF: ${RECOMBINATION_PA445}/${gene_name}.combined_region.einverted.genome_coords.gff \e[0m"


echo ""
echo "========================================================"
echo " PIPELINE COMPLETE - $(date)"
echo "========================================================"
echo "  Renamed FASTA   : ${filtered_assembly_renamed}"
echo "  CheckM          : ${CHECKM_PA445}/PA445.quality.checkm.tsv"
echo "  Bakta           : ${BAKTA_PA445}/"
echo "  Gene coordinates: ${GENE_PA445}/${gene_name}.coordinates.txt"
echo "  Gene sequences  : ${GENE_PA445}/${gene_name}_copy*.fasta"
echo "  Combined region : ${GENE_PA445}/${gene_name}.combined_region.fasta"
echo "  ISEScan region  : ${MGE_PA445}/${gene_name}.combined_region.isescan/"
echo "  ISEScan genome  : ${MGE_PA445}/PA445.full_genome.isescan/"
echo "  ISEScan GFF     : ${MGE_PA445}/${gene_name}.combined_region.isescan/${gene_name}.combined_region.isescan.genome_coords.unique.gff"
echo "  Integron region : ${MGE_PA445}/${gene_name}.combined_region.integron/"
echo "  Integron genome : ${MGE_PA445}/PA445.full_genome.integron/"
echo "  Integron GFF3   : ${MGE_PA445}/PA445.full_genome.integron/PA445.hybrid.filtered.integrons.gff3"
echo "  BLAST new_269   : ${BLAST_PA445}/new_269_vs_genome.txt"
echo "  Alignment       : ${MAFFT_PA445}/${gene_name}.aligned.fasta"
echo "  ICE_B81 BLAST           : ${EXT_PATH}/ice/ICE_B81_vs_PA445.txt"
echo "  ICEfinder submission    : ${EXT_PATH}/ice/PA445_for_ICEfinder_submission.fasta"
echo "  IS91 alignment          : ${EXT_PATH}/is91_boundaries/IS91_alignment.fasta"
echo "  IS91 self-BLAST         : ${EXT_PATH}/is91_boundaries/IS91_self_blast.txt"
echo "  Inverted repeats        : ${RECOMBINATION_PA445}/${gene_name}.combined_region.einverted.txt"
echo "  Palindromes             : ${RECOMBINATION_PA445}/${gene_name}.combined_region.palindrome.txt"
echo "  oriIS candidates        : ${RECOMBINATION_PA445}/PA445.combined_region.fuzznuc.IS91_oriIS_candidates.txt"
echo "  Synteny GBK inputs      : ${SYNTENY_PA445}/gbk_inputs/"
echo "  Synteny plot            : ${SYNTENY_PA445}/blaNDM1.combined_region.synteny.html"
echo "  PhiSpy results			: ${PHISPY_PA445}/PA445.full.prophages/"
echo "  Inverted repeat GFF     : ${RECOMBINATION_PA445}/${gene_name}.combined_region.einverted.genome_coords.gff"
echo "  Full log        : ${LOG}"
echo "========================================================"

echo " Manual steps still required:"
echo ""
echo -e "\e[32m 1. Alignment interpretation guide: \e[0m"
echo -e "\e[32m   Open ${MAFFT_PA445}/${gene_name}.aligned.fasta in AliView or MEGA \e[0m"
echo -e "\e[32m   100% identity = very recent IS91-mediated duplication \e[0m"
echo -e "\e[32m   <100% identity = older duplication with independent evolution \e[0m"
echo "  2. Update ICE_B81_ACCESSION at top of script with NCBI accession - OQ806931"
echo "     Need further confirmation"
echo "  3. Submit to ICEfinder: https://bioinfo-mml.sjtu.edu.cn/ICEfinder/"
echo "     Use file: ${EXT_PATH}/ice/PA445_for_ICEfinder_submission.fasta"
echo "  4. Add reference strain GBK files to synteny/gbk_inputs/ for clinker"
echo "     Confirmed blaNDM-1 carrying strains from RGI output are best candidates"
