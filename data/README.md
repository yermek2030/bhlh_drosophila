# data/

Processed inputs and intermediate results. Paths in the scripts are
relative to the repository root (`PROJECT_DIR`).

A few inputs are public datasets that are downloaded rather than copied here;
place them at the paths given:

| Path | Source |
|---|---|
| `chapter3/dm6.phyloP27way.bw` | UCSC, https://hgdownload.soe.ucsc.edu/goldenPath/dm6/phyloP27way/dm6.phyloP27way.bw |
| `chapter5/modERN_Resource_RNAi_Table_S43_trh_GENETICS-2022-305822.xlsx` | Fisher et al. (2023) *Genetics*, supplementary Table S43 (*trh* RNAi RNA-seq) |
| `chapter5/GSE28780_series_matrix.txt.gz`, `chapter5/GPL1322.soft.gz` | NCBI GEO, fetched by `GEOquery::getGEO()` in `chapter5_direct_targets.Rmd` |
| `vienna_tiles/vt_tiles_all_downloaded.csv` | Vienna Tiles (Kvon et al. 2014), https://enhancers.starklab.org/download/tiles/all ; see `vienna_tiles/PROVENANCE.txt` |
| `reference_sets/insitu_annot.csv` | BDGP in situ expression annotation, https://insitu.fruitfly.org/ |

Large intermediate files (genome-wide FASTA, shuffled sequence sets, genome-wide
FIMO scans, merged microarray tables) are written by the scripts named in the
README when they are run. Motif-analysis results are kept as the text outputs of
the MEME Suite and i-cisTarget (`meme.txt`, `streme.txt`, `tomtom.tsv`,
`fimo.tsv`, `statistics.tbl.gz`).

| Folder | Contents |
|---|---|
| `chapter1/` | Structure drawing (BioRender export) on which Figure 1.2 is labelled |
| `final_peaks/` | IDR 0.05 peak sets (narrowPeak) of Trh, Tgo, Sim, Dys, Sima, Twi, Vvl, Rib and Da; blacklist-filtered (`.clean`) sets; IDR consensus BED files and signal tables |
| `chapter3/` | Trh/Tgo joint peak set (614 peaks), peak FASTA (native and repeat-masked), Motif2 PWM `STRCACTGARARAAA.meme`, FIMO and STREME output, MEME-ChIP summaries, permutation and ChIP-Enrich results, histone-mark summaries |
| `chapter4/` | Sim/Tgo, Dys/Tgo and Sima/Tgo joint peak sets and FASTA, repeat-masked STREME runs, four-way TOMTOM run, FIMO output, 4×2 design tables, permutation results, CME-masking diagnostic |
| `chapter5/` | *trh* RNAi microarray tables (vsn/limma, R11 to R15, gzipped) and quality-control images, RT-qPCR data, direct-target tables, Ribbon and GSE28780 comparison tables |
| `ctrl_tf/` | FIMO output of Motif2 in the control-factor peak sets |
| `reference_gene_lists/` | EnhancerAtlas 2.0 embryonic enhancers (dm6), Vienna Tiles BED files and annotation, tissue reference gene sets, per-gene fold changes from Fisher et al. (2023) |
| `reference_sets/` | Tissue reference-set statistics (`refset_*.csv`) |
| `motif2_robustness/`, `motif2_gaga/`, `motif2_grh_confounds/`, `motif2_tagteam/`, `motif2_attribution/`, `motif2_concat_audit/` | Results of the Motif2 robustness, co-occurrence, attribution and concatenation-audit analyses |
| `icistarget/` | i-cisTarget region sets and result tables |
| `vienna_tiles/` | Source notes for the Vienna Tiles files |

Coordinates are dm6 with UCSC chromosome names. FIMO output on peak FASTA uses
peak-indexed coordinates (sequence name = peak identifier); genome-wide FIMO output
uses chromosomal coordinates.
