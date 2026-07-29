<h1>Whole genome analysis of <i>Pseudomonas aeruginosa</i> PA445 isolate and comparative genomic with global isolates</h1>

<h2>Background</h2>
	<p class="subtitle"><i>P. aeruginosa PA445</i> is a Vietnamese isolate carrying 2 blaNDM-1 gene copies, enhancing its resistance to beta-lactam drug class. Further annotations revealed surrounding flanking regions with mobile elements, proposing mechanism of action for its gene duplication.</p>
  <p class="subtitle">Comparative genomic of PA445 and other global isolates revealed similarities and distinction with other isolates carrying <i>blaNDM-1</i>, as well as suggestion of horizontal transfer of the <i>blaNDM-1</i> flanking region before duplication.</p>

<h2 id="conda-envs">Conda environments</h2>
	<table>
		<tr><th>Environment</th><th>Key tools</th></tr>
		<tr><td><code>preprocessing</code></td><td>Trimmomatic, FastQC, MultiQC</td></tr>
		<tr><td><code>assembly</code></td><td>SPAdes, QUAST, seqkit</td></tr>
		<tr><td><code>BPannotation</code></td><td>Bakta, Prokka, ABRicate</td></tr>
		<tr><td><code>BPstructure</code></td><td>CheckM</td></tr>
		<tr><td><code>BPtyping</code></td><td>mlst, pasty</td></tr>
		<tr><td><code>phylogeny</code></td><td>clinker, MAFFT, IQ-TREE3</td></tr>
		<tr><td><code>recombination</code></td><td>ISEScan, IntegronFinder, EMBOSS, PhiSpy</td></tr>
		<tr><td><code>pangenome</code></td><td>Panaroo</td></tr>
		<tr><td><code>ani</code></td><td>fastANI, ANIclustermap</td></tr>
    	<tr><td><code>ncbi</code></td><td>BLAST, NCBI datasets CLI, AMRFinderPlus, genbank_to</td></tr>
		<tr><td><code>rgi_env</code></td><td>RGI + CARD database</td></tr>
		<tr><td><code>plasmid</code></td><td>Platon, PlasmidFinder, MOB-suite</td></tr>
	</table>

<h2 id="repo-structure">Pipelines</h2>
	<table>
		<tr><th>Scripts</th><th>Purpose</th></tr>
		<tr><td><code>PA_annotation.bash</code></td><td>blaNDM-1 duplication-mechanism investigation</td></tr>
		<tr><td><code>PA_comparative.bash</code></td><td>Comparative genomics vs. PA445, genetic distance and time-caliberated phylogeny</td></tr>
	</table>

<div class="card">
		<div class="card-title">
			<h3><code>PA_annotation.bash</code></h3>
			<span class="tag">Bacterial &middot; single strain</span>
		</div>
		<div class="meta-row">
			<span><strong>Input:</strong> PA445 hybrid assembly (filtered, &ge;500&nbsp;bp contigs)</span>
			<span><strong>Focus gene:</strong> blaNDM-1</span>
		</div>
		<p>Annotates the PA445 genome with Bakta, extracts both blaNDM-1 copies and their flanking region, and investigates the duplication mechanism: insertion-sequence boundaries (ISEScan, IS91 alignment), integrons, the ICE_B81 mobile element, recombination/palindrome sites, prophage content (PhiSpy), and a clinker synteny figure.</p>
		<pre><code>./PA_annotation.bash</code></pre>
</div>

<div class="card">
		<div class="card-title">
			<h3><code>PA_comparative.bash</code></h3>
			<span class="tag">Bacterial &middot; comparative genomics</span>
		</div>
		<div class="meta-row">
			<span><strong>Input:</strong> PA445 + NCBI-downloaded global <em>P. aeruginosa</em> genomes</span>
			<span><strong>Stages:</strong> 8 phases</span>
		</div>
		<p>Downloads and filters global blaNDM-1 / background genomes, confirms true NDM-1 carriers (Prokka&nbsp;+&nbsp;RGI reconciliation), builds ST-aware subsampled and chromosomal-NDM-1-only datasets, and produces AMR heatmaps, core-gene phylogenies (Panaroo&nbsp;+&nbsp;IQ-TREE3), synteny plots, plasmid typing, and an optional BEAST-based time-calibrated phylogeny within PA445's lineage.</p>
		<pre><code>./PA_comparative.bash</code></pre>
</div>

<div class="card">
		<div class="card-title">
			<h3>Other python scripts</h3>
		</div>
		<div class="meta-row">
		<ul>
			<li><code>rename_contigs.py</code>: renaming contigs for consistencies before subsequent RGI profiling</li>
			<li><code>prokka_rgi_reconcile.py</code>: identify "true positive" calls for annotated antibiotic resistant genes if reported by both Prokka and RGI for an annotated gene</li>
			<li><code>gbk_region_extraction.py</code>: extraction of flanking region surrounding an annotated gene (e.g. <i>blaNDM-1</i>) with custom length (e.g. 50000bp)</li>
			<li><code>subsample_by_st.py</code>" subsampling dataset with conserved sequence type and inclusion of other background dataset optionally</li>
			<li><code>beast_prep_inputs_list.py</code> and <code>beast_make_xml.py</code>: preparation of inputs for BEAST time-caliberated phylogeny tree build</li>
		</ul>
</div>

<div class="card">
		<div class="card-title">
			<h3>Other Rscripts</h3>
		</div>
		<div class="meta-row">
		<ul>
			<li>phylogeny_O-type_AMR-heatmap folder: visualization of figure combining metadata + phylogeny tree + AMR profiling heatmap</li>
			<li><code>plot_time_tree.R</code>: visualization of time-caliberated phylogeny tree</li>
		</ul>
</div>

<h2 id="notes">Notes</h2>
	<ul>
		<li>Annotation format compatibility matters: Prokka GBK/GFF3 is the safe intermediate for downstream tools (Panaroo, IslandPath-DIMOB, GIPSy2); Bakta output frequently breaks them.</li>
		<li>SPAdes contig headers exceed some tools' LOCUS-name limits and are renamed to short sequential IDs before annotation.</li>
		<li>NCBI accession numbers referenced in scripts/comments should be independently verified against NCBI before citing.</li>
	</ul>

<footer>
		Internal lab pipelines &middot; run on <code>/storage/student9/</code> &middot; update this README when adding or restructuring a script.
</footer>
