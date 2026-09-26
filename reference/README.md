# reference/

Public reference files: download them into this folder (paths are relative to
the repository root):

| File | Source |
|---|---|
| `dm6.fa` | UCSC dm6 genome, https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz (gunzip) |
| `rmsk.txt`, `dm6_repeats.bed` | UCSC dm6 RepeatMasker table, https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz; `dm6_repeats.bed` is made from it by the first bash chunk of the motif section of `04_chapter_analyses/chapter3_trh_tgo_binding.Rmd` |
| `dm6_hardmasked.fa` | made from `dm6.fa` and `dm6_repeats.bed` with `bedtools maskfasta` (chunks in the Chapter 3 and Chapter 4 R Markdown files) |
| `liftover/dm2ToDm3.over.chain`, `liftover/dm3ToDm6.over.chain` | UCSC liftOver chains, https://hgdownload.soe.ucsc.edu/goldenPath/dm2/liftOver/ and https://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/ (gunzip) |
| `motif_databases/JASPAR/JASPAR2024_CORE_insects_non-redundant_v2.meme`, `motif_databases/JASPAR/JASPAR2024_CORE_insects_redundant_v2.meme`, `motif_databases/JASPAR/JASPAR2024_CORE_redundant_v2.meme`, `motif_databases/FLY/fly_factor_survey.meme`, `motif_databases/FLY/OnTheFly_2014_Drosophila.meme`, `motif_databases/FLY/flyreg.v2.meme`, `motif_databases/FLY/dmmpmm2009.meme`, `motif_databases/FLY/idmmpmm2009.meme` | MEME Suite motif database bundle, https://meme-suite.org/meme/db/motifs |
| `motif_databases/combined_drosophila.meme` | concatenation of `JASPAR2024_CORE_insects_non-redundant_v2.meme` and `fly_factor_survey.meme` (942 motifs), made by a bash chunk in `04_chapter_analyses/chapter3_trh_tgo_binding.Rmd` |
| `chipchip_histone_marks/*.bedGraph` | histone-mark ChIP-chip tracks (BDGP R4 / dm2 coordinates) used only by `04_chapter_analyses/chapter7/histone_marks_*.R` |

The ChIP-seq processing pipelines (`01_chipseq_processing/`) additionally need a
bowtie2 index of dm6 named `dm6` and the ENCODE blacklist
`dm6-blacklist.v2.bed` (https://github.com/Boyle-Lab/Blacklist) in their
working directory. The phyloP track is read from
`data/chapter3/dm6.phyloP27way.bw`
(https://hgdownload.soe.ucsc.edu/goldenPath/dm6/phyloP27way/dm6.phyloP27way.bw).
