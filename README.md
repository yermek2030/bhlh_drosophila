# Analysis code and processed data accompanying the PhD dissertation

> *Enhancer context and the shared binding architecture of bHLH-PAS-Tango
> heterodimers in late* Drosophila *embryogenesis.*

The dissertation analyses genome-wide binding of four heterodimers that share the
obligate partner Tango (Tgo): Trachealess/Tgo (Trh/Tgo), Single-minded/Tgo
(Sim/Tgo), Dysfusion/Tgo (Dys/Tgo) and Similar/Tgo (Sima/Tgo). The central result
is Motif2 (consensus `STRCACTGARARAAA`, 15 bp; PWM
`data/chapter3/STRCACTGARARAAA.meme`), an accessory motif enriched at distal
embryonic enhancers of all four heterodimers.

## Repository layout

| Folder | Contents |
|---|---|
| `01_chipseq_processing/` | ChIP-seq alignment and quality control, MACS3 threshold comparison, peak calling with IDR, pseudoreplicate IDR (Chapter 2) |
| `02_motif_analysis/` | Motif2 masking and rediscovery, checks of the STREME ranks and TOMTOM q-values, robustness tests, concatenation audit, co-occurrence (GAGA, Grainyhead) and i-cisTarget screens |
| `03_enrichment_and_permutation/` | Permutation tests (regioneR, ChIPpeakAnno) and ChIP-Enrich analyses; SLURM job files for the cluster runs |
| `04_chapter_analyses/` | R Markdown analyses for Chapters 3, 4 and 5, per-chapter scripts, and the figure scripts for Chapters 1, 3 and 4 |
| `05_structural_models/` | AlphaFold 3 models of the four heterodimers bound to a CME duplex, and the ChimeraX render script |
| `data/` | Processed inputs and intermediate results (peak sets, FASTA, FIMO/STREME/TOMTOM text output, tables); `data/README.md` lists the public inputs to download |
| `reference/` | Public reference files (genome, repeat annotation, motif databases), to be downloaded; see `reference/README.md` |
| `results/` | Figures, tables and `reported_values.tsv`, produced when the scripts are run |
| `environment/` | Software versions stated in the Methods chapter and the R packages used |
| `quick_checks/` | A short check of the joint peak sets against the counts in the text |

## Quick checks

From the repository root (R 4 with GenomicRanges, GenomeInfoDb and rtracklayer):

```sh
export PROJECT_DIR=$(pwd)
Rscript quick_checks/joint_peak_count.R           # Trh/Tgo 614; Sim/Tgo 355; Dys/Tgo 870; Sima/Tgo 296
Rscript 02_motif_analysis/verify_tomtom_worst_q.R  # four-way TOMTOM comparison of the Motif2 PWMs
```

`verify_tomtom_worst_q.R` identifies, among the four-way TOMTOM runs in `data/chapter4/`, the
one described in Chapter 4 (`data/chapter4/tomtom_motif2_4way/`; worst pair
Sima then Sim, best pair Dys then Trh) and reports its worst q-value over the 12
cross-comparisons (2.6e-10, the value given in the text).

```sh
bash 02_motif_analysis/verify_streme_ranks.sh      # STREME rank, E-value and consensus of Motif2 per heterodimer
```

## Running the analyses

* Set `PROJECT_DIR` to the repository root, or run R from it. Every script and
  R Markdown file is self-contained: it defines the paths, colours and small
  functions it uses at the top, and reads only the files under `data/`. All
  coordinates are dm6 (UCSC chromosome names); analyses are restricted to
  `chr2L chr2R chr3L chr3R chr4 chrX`.
* The R Markdown files are written to be run chunk by chunk, in order, in one R
  session (chunks carry `eval=FALSE`, so knitting does not execute them). Standalone
  `.R` scripts are run with `Rscript` from the repository root.
* Long permutation and ChIP-Enrich runs were run on a SLURM cluster: the `.slurm`
  files hold the resources used; set the partition and account for your cluster.
* Each script writes the values it contributes to the text to
  `results/reported_values.tsv` (columns `section`, `name`, `value`); a value that
  is already there is replaced, so scripts can be re-run singly and in any order.
* Motif analyses use MEME Suite 5.5.5; `meme-chip`, `streme`, `fimo`, `tomtom` and
  `fasta-shuffle-letters` must be on `PATH` (`ch4_motif2_attribution.R` also
  accepts `MEME_BIN`, the MEME `bin/` directory).
* The ChIP-seq processing pipelines (`01_chipseq_processing/`) run in a separate
  working directory that holds the FASTQ files and receives `03_bam/` to `11_final/`;
  their final peak sets (`11_final/*_IDR0.05.final.narrowPeak`) and the MACS3
  summit files of the Trh replicates are in `data/final_peaks/`. They
  need a bowtie2 index of dm6 named `dm6` and `dm6-blacklist.v2.bed` in that
  directory. The deepTools TSS-profile chunk of `chapter3_trh_tgo_binding.Rmd`
  also runs there, on the blacklist-filtered BAM files.

Suggested order: `chapter3_trh_tgo_binding.Rmd`, then `chapter3/ch3_peak_annotation_counts.R`
and the other `chapter3/` scripts; `chapter4_heterodimer_comparison.Rmd`, then
`chapter4_control_factors.Rmd`, `chapter4/ch4_ctrl_motif2_panel.R` and the other
`chapter4/` scripts (`ch4_motif2_attribution.R` and `ch4_ctrl_motif2_symmetric.R`
before `02_motif_analysis/motif2_dm6_background.R`); `chapter5_direct_targets.Rmd`
(self-contained), then `chapter5/ch5_rule_on_control_factors.R`; finally
`chapter3_4_figures.Rmd` (figures; reads stored results only). In
`02_motif_analysis/robustness/`, run `motif2_robustness_b_heldout.R` before
`motif2_robustness_d_heldout_halfA_pwm.R`; both it and `motif2_dm6_background.R`
read the shuffled sequences that `motif2_robustness_b_heldout.R` writes.

## Where each part of the dissertation is computed

| Dissertation | Code | Main outputs |
|---|---|---|
| Ch 2, ChIP-seq alignment and QC (bowtie2, Picard, blacklist filtering, SPP cross-correlation, NSC/RSC) | `01_chipseq_processing/alignment_qc/alignment_qc_{trh,tgo,sim,dys,sima}.Rmd`, `alignment_qc_twi_kr_stat92e.Rmd`, `alignment_qc_vvl_rib.Rmd`, `spp_cross_correlation_vvl_rib.sh`; for Da the alignment, QC and peak calling are one script, `peak_calling_idr/production_peaks_da.sh` | BAM files, QC tables (processing directory) |
| Ch 2, MACS3 p-value threshold comparison | `01_chipseq_processing/macs3_threshold_comparison/threshold_comparison_*.sh` | peak counts per threshold |
| Ch 2, peak calling, IDR and consensus peak sets | `01_chipseq_processing/peak_calling_idr/production_peaks_*.sh` (one per factor; Twi, Vvl/Rib and Da are the controls), `final_peak_set_{trh,tgo,sim}.sh` | `data/final_peaks/*_IDR0.05.final.narrowPeak` |
| Ch 2, file-level facts about the Motif2 inputs (MEME-ChIP centring, masked bases, PWM header and origin, FIMO background) | `02_motif_analysis/verify_motif2_provenance.R` | `results/reported_values.tsv` |
| Ch 2, pseudoreplicate IDR (Trh, Sim, Sima) | `01_chipseq_processing/pseudoreplicates/pseudoreplicate_idr_*.sh` | pseudoreplicate IDR counts |
| Ch 2, microarray quality control (Dapple spot classes, dye-swap orientation) | `04_chapter_analyses/chapter2/dapple_spot_class_qc.R`, `dye_swap_orientation_check.R` | `results/tables/dapple_spot_class_by_stage.tsv` |
| Ch 3, Trh/Tgo binding: joint peak set (614), genomic distribution, TSS profiles, annotation, overlap statistics, enhancer and promoter context, reference gene sets, ChIP-Enrich, CME stratification, MEME-ChIP triage, FIMO, Motif2 enrichment, regression and conservation tests | `04_chapter_analyses/chapter3_trh_tgo_binding.Rmd` | `data/chapter3/`, `results/` |
| Ch 3, occupancy on the input-normalised scale: promoter against distal-enhancer fold-enrichment, and the regression of fold-enrichment on Motif2 dosage (writes the fitted models read by `ch3_tab_motif2_regression.R`) | `04_chapter_analyses/chapter3/ch3_occupancy_fold_enrichment.R` | `data/chapter3/test1_regression_binding.rds`, `results/reported_values.tsv` |
| Ch 3, counts and statistics computed from the primary files | `04_chapter_analyses/chapter3/ch3_peak_annotation_counts.R`, `ch3_memechip_inventory.R`, `ch3_tab_motif2_regression.R`, `ch3_tab_motif2_gaga.R`, `ch3_motif2_auroc.R`, `ch3_motif2_summit_annulus_stats.R`, `ch3_motif2_tagteam.R`, `ch3_threshold_sensitivity_values.R`, `ch3_fig_cme_stratification.R` | `results/reported_values.tsv`, `results/tables/`, `results/figures/` |
| Ch 3 and 4, tissue reference gene set: composition, sharing, pairwise, tissue-contrast and power statistics | `04_chapter_analyses/chapter3/ch3_ch4_reference_set_statistics.R` | `data/reference_sets/refset_*.csv` |
| Ch 3, gene-length-matched null for the tracheal/salivary-gland reference set and its partitions | `04_chapter_analyses/chapter3/ch3_tissue_enrichment_length_matched.R` | `data/reference_sets/refset_lenmatched_5000.csv` |
| Ch 3 and 4, composition nulls for CME carriage (dinucleotide shuffles, random intervals) | `04_chapter_analyses/chapter4/ch4_cme_composition_null.R` | `data/chapter4/cme_composition_null/` |
| Ch 3, permutation tests (Trh×Tgo overlap, TSS-matched null, enhancer overlap, tissue-TSS proximity, Motif2-enhancer null) | `03_enrichment_and_permutation/permutation_tests/perm_overlap_trh_tgo.R`, `perm_overlap_trh_tgo_tss_matched.R`, `perm_enhancer_overlap.R`, `perm_enhancers_tss_matched_trh_tgo.R`, `perm_tissue_tss_proximity.R`, `perm_motif2_enhancer_null.R` | `data/chapter3/perm_*.rds`, `*.csv` |
| Ch 3, ChIP-Enrich (tissue gene set, GO; proximal/distal split) | `03_enrichment_and_permutation/chipenrich/chipenrich_trh_tgo.R`, `chipenrich_trh_tgo_distal_split.R` | `data/chapter3/chipenrich_ch3*.rds`, `*.csv` |
| Ch 3, Motif2 robustness (shuffled null, held-out halves, register, summit centrality) | `02_motif_analysis/robustness/` | `data/motif2_robustness/` |
| Ch 3 and 4, held-out half scanned with the matrix learned on the other half | `02_motif_analysis/robustness/motif2_robustness_d_heldout_halfA_pwm.R` | `data/motif2_robustness/heldout_halfA_pwm/` |
| Ch 3, Motif2 co-occurrence with GAGA and Grainyhead sites; i-cisTarget screens | `02_motif_analysis/cooccurrence/`, `02_motif_analysis/icistarget/` | `data/motif2_gaga/`, `data/motif2_grh_confounds/`, `data/icistarget/` |
| Ch 2 and 4, CME-masked rediscovery of Motif2 | `02_motif_analysis/cme_masking_rediscovery.sh` | `data/chapter4/motif2_diagnostic/` |
| Ch 4, Sim/Tgo, Dys/Tgo and Sima/Tgo joint peak sets, enhancer overlap, 4×2 tissue design, regulatory classes, CME, Motif2 across heterodimers, TOMTOM equivalence | `04_chapter_analyses/chapter4_heterodimer_comparison.Rmd` | `data/chapter4/`, `results/tables/ch4_motif2_tomtom_pairwise.tsv` |
| Ch 4 and appendix, repeat-masked STREME per heterodimer and rank listing | `04_chapter_analyses/chapter4_heterodimer_comparison.Rmd` (STREME chunk), `02_motif_analysis/verify_streme_ranks.sh` | `data/chapter3/Trh_Tgo_joint_peaks_filt_repeatmasked_streme/`, `data/chapter4/*_joint_peaks_filt_repeatmasked_streme/` |
| Ch 4, permutation tests and ChIP-Enrich for the four heterodimers | `03_enrichment_and_permutation/permutation_tests/perm_overlap_heterodimers.R`, `perm_4x2.R`, `perm_enhancers_tss_matched_all.R`; `03_enrichment_and_permutation/chipenrich/chipenrich_4x2.R` | `data/chapter4/table4_*.csv`, `perm_*.rds` |
| Ch 4, control factors (Twi, Vvl, Rib, Da) and the Motif2 control panel | `04_chapter_analyses/chapter4_control_factors.Rmd`, `chapter4/ch4_ctrl_motif2_panel.R`, `ch4_da_motif2.R`, `ch4_vvl_trh_cooccupancy_test.R`, `ch4_singlefactor_motif2.R` | `data/ctrl_tf/`, `data/chapter4/table4_ctrl_motif2_frequency.csv` |
| Ch 4, Motif2 attribution (bound vs unbound enhancers, non-enhancer baseline) | `04_chapter_analyses/chapter4/ch4_motif2_attribution.R`, `ch4_motif2_nonenh_baseline.R` | `data/motif2_attribution/` |
| Ch 4, heterodimer × tissue specificity test | `04_chapter_analyses/chapter4/ch4_specificity_test.R` | `data/chapter4/table4_specificity_disjoint.csv` |
| Ch 4, gene-length-matched null for the eight cells of the 4×2 matrix | `04_chapter_analyses/chapter4/ch4_tissue_enrichment_length_matched_4x2.R` | `data/chapter4/lenmatch_4x2/` |
| Ch 4, control factors processed in the same way as the Trh/Tgo reference (alone, with Tgo, without Tgo; window-GC-matched background) | `04_chapter_analyses/chapter4/ch4_ctrl_motif2_symmetric.R` | `data/ctrl_tf/symmetric/` |
| Ch 4, de novo discovery on peaks exclusive to one heterodimer | `04_chapter_analyses/chapter4/ch4_exclusive_peak_discovery.R` | `data/chapter4/exclusive_peak_discovery/` |
| Ch 4, CME occurrence per heterodimer | `04_chapter_analyses/chapter4/ch4_cme_occurrence.R` | `results/reported_values.tsv` |
| Ch 4, Vienna Tiles validation | `03_enrichment_and_permutation/permutation_tests/prepare_vienna_tiles_dm6.R`, `perm_enhancers_vienna_tiles.R`; `04_chapter_analyses/chapter4/ch4_vienna_tiles_summary.R` | `data/reference_gene_lists/vt_*.bed`, `data/chapter4/` |
| Ch 4 and appendix, concatenation audit of Motif2 (information content, background-order sweep and its E-values, register, core-only tests) | `02_motif_analysis/concat_audit/` (run `RUN_ALL.sh`), `04_chapter_analyses/chapter4/ch4_fig_concat_audit.R` | `data/motif2_concat_audit/` |
| Appendix, Motif2 prevalence comparisons rescanned against a dm6 background | `02_motif_analysis/motif2_dm6_background.R` | `data/motif2_dm6_background/` |
| Ch 5, direct targets (ChIP-seq intersected with trh RNAi RNA-seq), sensitivity to the symbol mapping, microarray concordance, RT-qPCR, Ribbon and GSE28780 comparisons | `04_chapter_analyses/chapter5_direct_targets.Rmd` | `data/chapter5/`, `results/figures/` |
| Ch 5, microarray: probes passing Benjamini-Hochberg thresholds per stage and the minimum detectable effect at stage 13 | `04_chapter_analyses/chapter5/ch5_microarray_stage_counts.R` | `results/reported_values.tsv` |
| Ch 5, the direct-target rule applied to the control factors and the other heterodimers | `04_chapter_analyses/chapter5/ch5_rule_on_control_factors.R` | `data/chapter5/ctrl_rule_benchmark/` |
| Ch 3 and 4 figures | `04_chapter_analyses/chapter3_4_figures.Rmd` (last chunk rebuilds all), `build_chapter3_4_figures.R` | `results/figures/` |
| Ch 7 and appendix, histone-mark ChIP-chip analysis | `04_chapter_analyses/chapter7/histone_marks_{1..4}_*.R` | `data/chapter3/histone_marks_*` |
| Ch 1, Figure 1.2 (heterodimer schematic: labels, leader lines and CME pentamer strip drawn on the BioRender export) | `04_chapter_analyses/chapter1/fig1_2_bhlh_pas_architecture.R` | `results/figures/fig1_2_bhlh_pas_architecture.png` |
| Ch 1 (methods in Ch 2), structural models of the heterodimers (AlphaFold 3; seeds Trh 456, Sim 999, Sima 42, Dys 777) | `05_structural_models/render_af3_models.cxc` | `05_structural_models/models/`, `results/figures/panel_*_aligned.png` |

## Input datasets

ChIP-seq data are from the ENCODE portal (modERN), embryo stage as described in
Chapter 2:

| Factor | Experiment | Replicate files | Control experiment | Control files |
|---|---|---|---|---|
| Trh | ENCSR153KTS | ENCFF334EBW, ENCFF479LGT | ENCSR235TIB | ENCFF022NFX |
| Tgo | ENCSR306TVS | ENCFF726TZP, ENCFF500UCX, ENCFF124RRH, ENCFF516EGV, ENCFF724VYL, ENCFF270GLQ | ENCSR400MGD | ENCFF098LKV, ENCFF305YAC |
| Sim | ENCSR436KCT | ENCFF946PGU, ENCFF965DLV, ENCFF130YRG, ENCFF586UQK | ENCSR449OOC | ENCFF987WFH, ENCFF098IAX |
| Sima | ENCSR390QUD | ENCFF048IYN, ENCFF309ZFF | ENCSR802IDG | ENCFF936AOH |
| Dys | ENCSR522KTT | ENCFF850GFV, ENCFF038CZD, ENCFF043NYJ, ENCFF362XNT, ENCFF207TJL, ENCFF027UZF | ENCSR361NFK | ENCFF063EMW, ENCFF833TRM |
| Twi | ENCSR119OUI | ENCFF918UHG, ENCFF271FJU | ENCSR981YAA | ENCFF082BHL |
| Vvl | ENCSR350VLI | ENCFF895FAF, ENCFF061PMX, ENCFF332TNC, ENCFF626ISZ | ENCSR186OZX | ENCFF218OYW, ENCFF466TTY |
| Rib | ENCSR129UZN | ENCFF596NQN, ENCFF929RDR, ENCFF489JKI, ENCFF889PVP | ENCSR600RUY | ENCFF458LQH, ENCFF037LQP |
| Da | ENCSR416MJQ | ENCFF303EVV, ENCFF831FKP | ENCSR515BLI | ENCFF405BXW |

Candidate controls evaluated and excluded (Chapter 4): Krüppel (ENCSR156QJC),
Stat92E (ENCSR240ADR, ENCSR290OJD, ENCSR494UZY) and Tramtrack (ENCSR362ZUH).

Other inputs: *trh* RNAi RNA-seq of Fisher et al. (2023), Table S43 (per-gene fold
changes in `data/reference_gene_lists/Fisher2023_per_gene_fold_changes.csv`; the
original workbook is the journal's supplementary Table S43); the in-house *trh* RNAi FlyChip
microarray (vsn/limma tables `data/chapter5/R11` to `R15.vsn.limma.fdr.txt.gz`, RT-qPCR
in `data/chapter5/`); GEO GSE28780 (Chung et al. 2011, platform GPL1322); embryonic
enhancers of EnhancerAtlas 2.0 (Gao and Qian 2020;
`data/reference_gene_lists/merged_embryo_5-20_enhancers_dm6.bed`); Vienna Tiles
(Kvon et al. 2014); Ribbon direct targets (Loganathan et al. 2016;
`data/chapter5/rib-ts-sg-direct-targets-complete.csv`); BDGP in situ annotation;
dm6 genome, RepeatMasker, phyloP27way and liftOver chains (UCSC); JASPAR 2024 insect
and Fly Factor Survey motif collections. Download locations are listed in
`reference/README.md` and `data/README.md`.

## Software

Tool versions are those stated in the Methods chapter
(`environment/software_versions.tsv`); `environment/r_packages.tsv` lists every R
package used. Cluster jobs used the conda environments `phd` (general),
`idr` (IDR 2.0.4.2) and `r_phylo` (SPP); motif analyses used MEME Suite 5.5.5.

