#!/usr/bin/env python3
"""
make_beast_xml.py
=================
Generate a BEAST 2.7 tip-dated XML from a (recombination-free) SNP alignment and
a decimal-date table, so the lineage can be dated non-interactively.

Model (standard within-lineage bacterial setup, all overridable in BEAUti):
  * GTR + Gamma(4) site model
  * relaxed lognormal molecular clock (UCLD)   [--clock strict for a strict clock]
  * Coalescent Constant Population tree prior   [--tree-prior skyline for BSP]
  * tip dates forward in time (most recent = largest value)

IMPORTANT
---------
This XML is a STARTING point. It is emitted as well-formed BEAST 2.7 XML, but it
is not validated against a running BEAST here. Before a long run:
  1. open it in BEAUti (File > Load) or run a short smoke chain to confirm it
     loads in *your* BEAST version (2.6 vs 2.7 namespaces differ);
  2. if your data are SNP-only, set the ascertainment / constant-site correction
     (the constant A,C,G,T counts are written to a sidecar by the pipeline);
  3. after running, check ESS > 200 for all parameters in Tracer before trusting
     any date - convergence is not automatic.

Usage:
  make_beast_xml.py --alignment snps.fasta --dates dates.tsv --out beast.xml \
      [--clock relaxed|strict] [--tree-prior constant|skyline] \
      [--chain-length 100000000] [--log-every 10000] [--prefix run1]
"""

import argparse
import os
import sys
from xml.sax.saxutils import escape


def eprint(*a, **k):
    print(*a, file=sys.stderr, **k)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def read_dates(path):
    dates = {}
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h.lower(): i for i, h in enumerate(header)}
        ti, di = idx.get("taxon", 0), idx.get("decimal_date", 1)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > max(ti, di):
                dates[cols[ti]] = cols[di]
    return dates


def parse_args(argv=None):
    ap = argparse.ArgumentParser(prog="make_beast_xml.py",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description="Generate a BEAST 2.7 tip-dated XML.")
    ap.add_argument("--alignment", required=True, help="SNP FASTA (recombination-free)")
    ap.add_argument("--dates", required=True, help="TSV with taxon, decimal_date columns")
    ap.add_argument("--out", required=True, help="output BEAST XML")
    ap.add_argument("--clock", choices=["relaxed", "strict"], default="relaxed")
    ap.add_argument("--tree-prior", choices=["constant", "skyline"], default="constant")
    ap.add_argument("--chain-length", type=int, default=100_000_000)
    ap.add_argument("--log-every", type=int, default=10_000)
    ap.add_argument("--prefix", default=None, help="log/tree file prefix (default: out basename)")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if not os.path.isfile(args.alignment):
        eprint(f"ERROR: alignment not found: {args.alignment}"); sys.exit(2)
    if not os.path.isfile(args.dates):
        eprint(f"ERROR: dates not found: {args.dates}"); sys.exit(2)

    seqs = read_fasta(args.alignment)
    dates = read_dates(args.dates)
    taxa = [t for t in seqs if t in dates]
    missing = [t for t in seqs if t not in dates]
    if missing:
        eprint(f"WARNING: {len(missing)} aligned taxa have no date and are dropped: "
               f"{', '.join(missing[:5])}{'...' if len(missing) > 5 else ''}")
    if len(taxa) < 3:
        eprint(f"ERROR: only {len(taxa)} dated taxa in the alignment; too few to date."); sys.exit(2)

    prefix = args.prefix or os.path.splitext(os.path.basename(args.out))[0]
    seq_blocks = "\n".join(
        f'        <sequence id="seq_{escape(t)}" spec="Sequence" taxon="{escape(t)}" '
        f'totalcount="4" value="{seqs[t]}"/>' for t in taxa)
    trait = ",\n".join(f"{t}={dates[t]}" for t in taxa)

    # ---- clock block ----
    if args.clock == "relaxed":
        clock_state = (
            '        <parameter id="ucldMean.c" spec="parameter.RealParameter" name="stateNode">1.0E-6</parameter>\n'
            '        <parameter id="ucldStdev.c" spec="parameter.RealParameter" lower="0.0" name="stateNode">0.1</parameter>\n'
            '        <stateNode id="rateCategories.c" spec="parameter.IntegerParameter" dimension="10">1</stateNode>')
        branch_rate = (
            '            <branchRateModel id="RelaxedClock.c" spec="beast.base.evolution.branchratemodel.UCRelaxedClockModel" '
            'clock.rate="@ucldMean.c" rateCategories="@rateCategories.c" tree="@Tree.t">\n'
            '                <LogNormal id="LogNormalDistributionModel.c" S="@ucldStdev.c" meanInRealSpace="true" name="distr">\n'
            '                    <parameter id="RealParameter.M" spec="parameter.RealParameter" estimate="false" name="M">1.0</parameter>\n'
            '                </LogNormal>\n'
            '            </branchRateModel>')
        clock_priors = (
            '            <prior id="MeanRatePrior.c" name="distribution" x="@ucldMean.c">\n'
            '                <LogNormal id="LogNormal.0" name="distr" M="-13.8" S="2.0"/>\n'
            '            </prior>\n'
            '            <prior id="ucldStdevPrior.c" name="distribution" x="@ucldStdev.c">\n'
            '                <Gamma id="Gamma.0" name="distr" alpha="0.5396" beta="0.3819"/>\n'
            '            </prior>')
        clock_ops = (
            '        <operator id="ucldMeanScaler.c" spec="kernel.BactrianScaleOperator" parameter="@ucldMean.c" weight="1.0"/>\n'
            '        <operator id="ucldStdevScaler.c" spec="kernel.BactrianScaleOperator" parameter="@ucldStdev.c" weight="3.0"/>\n'
            '        <operator id="CategoriesRandomWalk.c" spec="operator.IntRandomWalkOperator" parameter="@rateCategories.c" weight="10.0" windowSize="1"/>\n'
            '        <operator id="CategoriesSwapOperator.c" spec="operator.SwapOperator" intparameter="@rateCategories.c" weight="10.0"/>\n'
            '        <operator id="CategoriesUniform.c" spec="operator.UniformOperator" parameter="@rateCategories.c" weight="10.0"/>')
        clock_log = '        <log idref="ucldMean.c"/>\n        <log idref="ucldStdev.c"/>'
    else:
        clock_state = '        <parameter id="clockRate.c" spec="parameter.RealParameter" name="stateNode">1.0E-6</parameter>'
        branch_rate = ('            <branchRateModel id="StrictClock.c" spec="beast.base.evolution.branchratemodel.StrictClockModel" '
                       'clock.rate="@clockRate.c"/>')
        clock_priors = ('            <prior id="ClockPrior.c" name="distribution" x="@clockRate.c">\n'
                        '                <LogNormal id="LogNormal.0" name="distr" M="-13.8" S="2.0"/>\n'
                        '            </prior>')
        clock_ops = '        <operator id="clockScaler.c" spec="kernel.BactrianScaleOperator" parameter="@clockRate.c" weight="3.0"/>'
        clock_log = '        <log idref="clockRate.c"/>'

    # ---- tree prior ----
    if args.tree_prior == "skyline":
        tp_state = ('        <parameter id="bPopSizes.t" spec="parameter.RealParameter" dimension="5" lower="0.0" name="stateNode">1.0</parameter>\n'
                    '        <stateNode id="bGroupSizes.t" spec="parameter.IntegerParameter" dimension="5">1</stateNode>')
        tp_dist = ('            <distribution id="BayesianSkyline.t" spec="beast.base.evolution.tree.coalescent.BayesianSkyline" '
                   'groupSizes="@bGroupSizes.t" popSizes="@bPopSizes.t" tree="@Tree.t"/>')
        tp_priors = ('            <prior id="popSizesPrior.t" name="distribution" x="@bPopSizes.t">\n'
                     '                <OneOnX id="OneOnX.0" name="distr"/>\n'
                     '            </prior>')
        tp_ops = ('        <operator id="popSizesScaler.t" spec="kernel.BactrianScaleOperator" parameter="@bPopSizes.t" weight="15.0"/>\n'
                  '        <operator id="groupSizesDelta.t" spec="operator.DeltaExchangeOperator" integer="true" weight="6.0">\n'
                  '            <intparameter idref="bGroupSizes.t"/>\n'
                  '        </operator>')
        tp_log = '        <log idref="BayesianSkyline.t"/>\n        <log idref="bPopSizes.t"/>'
    else:
        tp_state = '        <parameter id="popSize.t" spec="parameter.RealParameter" name="stateNode">1.0</parameter>'
        tp_dist = ('            <distribution id="CoalescentConstant.t" spec="Coalescent">\n'
                   '                <populationModel id="ConstantPopulation.t" spec="ConstantPopulation" popSize="@popSize.t"/>\n'
                   '                <treeIntervals id="TreeIntervals.t" spec="beast.base.evolution.tree.TreeIntervals" tree="@Tree.t"/>\n'
                   '            </distribution>')
        tp_priors = ('            <prior id="PopSizePrior.t" name="distribution" x="@popSize.t">\n'
                     '                <OneOnX id="OneOnX.0" name="distr"/>\n'
                     '            </prior>')
        tp_ops = '        <operator id="PopSizeScaler.t" spec="kernel.BactrianScaleOperator" parameter="@popSize.t" weight="3.0"/>'
        tp_log = '        <log idref="popSize.t"/>'

    xml = f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast beautitemplate='Standard' beautistatus='' version="2.7"
       namespace="beast.base.evolution.alignment:beast.base.core:beast.base.evolution.tree.coalescent:beast.base.core.util:beast.base.evolution.nuc:beast.base.evolution.operator:beast.base.evolution.operator.kernel:beast.base.inference.operator:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood:beast.base.evolution.branchratemodel:beast.base.inference:beast.base.inference.parameter:beast.base.evolution.tree:beast.base.inference.distribution">

    <data id="alignment" spec="Alignment" dataType="nucleotide">
{seq_blocks}
    </data>

    <map name="Uniform">beast.base.inference.distribution.Uniform</map>
    <map name="LogNormal">beast.base.inference.distribution.LogNormalDistributionModel</map>
    <map name="Gamma">beast.base.inference.distribution.Gamma</map>
    <map name="OneOnX">beast.base.inference.distribution.OneOnX</map>
    <map name="prior">beast.base.inference.distribution.Prior</map>
    <map name="Exponential">beast.base.inference.distribution.Exponential</map>

    <run id="mcmc" spec="MCMC" chainLength="{args.chain_length}">
        <state id="state" spec="State" storeEvery="5000">
            <tree id="Tree.t" spec="beast.base.evolution.tree.Tree" name="stateNode">
                <trait id="dateTrait" spec="beast.base.evolution.tree.TraitSet" traitname="date-forward" value="{trait}">
                    <taxa id="TaxonSet" spec="TaxonSet" alignment="@alignment"/>
                </trait>
                <taxonset idref="TaxonSet"/>
            </tree>
            <parameter id="gammaShape.s" spec="parameter.RealParameter" name="stateNode">1.0</parameter>
            <parameter id="rateAC.s" spec="parameter.RealParameter" name="stateNode">1.0</parameter>
            <parameter id="rateAG.s" spec="parameter.RealParameter" name="stateNode">1.0</parameter>
            <parameter id="rateAT.s" spec="parameter.RealParameter" name="stateNode">1.0</parameter>
            <parameter id="rateCG.s" spec="parameter.RealParameter" name="stateNode">1.0</parameter>
            <parameter id="rateGT.s" spec="parameter.RealParameter" name="stateNode">1.0</parameter>
            <parameter id="freqParameter.s" spec="parameter.RealParameter" dimension="4" lower="0.0" upper="1.0" name="stateNode">0.25</parameter>
{clock_state}
{tp_state}
        </state>

        <init id="RandomTree.t" spec="RandomTree" estimate="false" initial="@Tree.t" taxa="@alignment">
            <populationModel id="InitPop" spec="ConstantPopulation">
                <parameter id="initPopSize" spec="parameter.RealParameter" name="popSize">1.0</parameter>
            </populationModel>
        </init>

        <distribution id="posterior" spec="CompoundDistribution">
            <distribution id="prior" spec="CompoundDistribution">
{tp_dist}
                <prior id="GammaShapePrior.s" name="distribution" x="@gammaShape.s">
                    <Exponential id="Exponential.0" name="distr" mean="1.0"/>
                </prior>
{clock_priors}
{tp_priors}
            </distribution>
            <distribution id="likelihood" spec="CompoundDistribution" useThreads="true">
                <distribution id="treeLikelihood" spec="TreeLikelihood" data="@alignment" tree="@Tree.t">
                    <siteModel id="SiteModel.s" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape.s">
                        <substModel id="gtr" spec="GTR" rateAC="@rateAC.s" rateAG="@rateAG.s" rateAT="@rateAT.s" rateCG="@rateCG.s" rateGT="@rateGT.s">
                            <frequencies id="estimatedFreqs.s" spec="Frequencies" frequencies="@freqParameter.s"/>
                        </substModel>
                    </siteModel>
{branch_rate}
                </distribution>
            </distribution>
        </distribution>

        <operator id="treeScaler.t" spec="kernel.BactrianScaleOperator" scaleFactor="0.5" tree="@Tree.t" weight="3.0"/>
        <operator id="treeRootScaler.t" spec="kernel.BactrianScaleOperator" rootOnly="true" scaleFactor="0.7" tree="@Tree.t" weight="3.0"/>
        <operator id="UniformOperator.t" spec="kernel.BactrianNodeOperator" tree="@Tree.t" weight="30.0"/>
        <operator id="SubtreeSlide.t" spec="kernel.BactrianSubtreeSlide" tree="@Tree.t" weight="15.0"/>
        <operator id="narrow.t" spec="Exchange" tree="@Tree.t" weight="15.0"/>
        <operator id="wide.t" spec="Exchange" isNarrow="false" tree="@Tree.t" weight="3.0"/>
        <operator id="WilsonBalding.t" spec="WilsonBalding" tree="@Tree.t" weight="3.0"/>
        <operator id="gammaShapeScaler.s" spec="kernel.BactrianScaleOperator" parameter="@gammaShape.s" weight="0.1"/>
        <operator id="RateACScaler.s" spec="kernel.BactrianScaleOperator" parameter="@rateAC.s" weight="0.1"/>
        <operator id="RateAGScaler.s" spec="kernel.BactrianScaleOperator" parameter="@rateAG.s" weight="0.1"/>
        <operator id="RateATScaler.s" spec="kernel.BactrianScaleOperator" parameter="@rateAT.s" weight="0.1"/>
        <operator id="RateCGScaler.s" spec="kernel.BactrianScaleOperator" parameter="@rateCG.s" weight="0.1"/>
        <operator id="RateGTScaler.s" spec="kernel.BactrianScaleOperator" parameter="@rateGT.s" weight="0.1"/>
        <operator id="FrequenciesExchanger.s" spec="kernel.BactrianDeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s"/>
{clock_ops}
{tp_ops}

        <logger id="tracelog" spec="Logger" fileName="{prefix}.log" logEvery="{args.log_every}" model="@posterior">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log idref="treeLikelihood"/>
            <log id="TreeHeight.t" spec="beast.base.evolution.tree.TreeStatLogger" tree="@Tree.t"/>
            <log idref="gammaShape.s"/>
{clock_log}
{tp_log}
        </logger>
        <logger id="screenlog" spec="Logger" logEvery="{args.log_every}">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
        </logger>
        <logger id="treelog.t" spec="Logger" fileName="{prefix}.trees" logEvery="{args.log_every}" mode="tree">
            <log id="TreeWithMetaDataLogger.t" spec="beast.base.evolution.TreeWithMetaDataLogger" tree="@Tree.t"/>
        </logger>
    </run>
</beast>
'''

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(xml)

    eprint("==================== make_beast_xml ====================")
    eprint(f"  taxa in XML   : {len(taxa)}")
    eprint(f"  clock         : {args.clock}")
    eprint(f"  tree prior    : {args.tree_prior}")
    eprint(f"  chain length  : {args.chain_length:,}")
    eprint(f"  -> {args.out}  (log/trees prefix: {prefix})")
    eprint("  REMINDER: verify it loads in BEAUti and check ESS>200 in Tracer.")
    eprint("=======================================================")
    return 0


if __name__ == "__main__":
    sys.exit(main())
