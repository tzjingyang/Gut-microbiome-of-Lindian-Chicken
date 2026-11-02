# Analysis Pipeline Summary

**Study**: Muscle fiber diameter (MFD) in Lindian chickens — host whole-genome GWAS integrated with 16S rRNA microbiome profiling across five intestinal segments (duodenum, jejunum, ileum, cecum, rectum; 300 birds × 5 segments).

All analyses were run in R 4.3.3 within two conda environments: a host-GWAS/statistics environment and a microbiome-plotting environment (see version table below).

---

## 1. Software and reference databases

| Software / database | Version | Used for |
|---|---|---|
| R | 4.3.3 | All statistical analyses and figure rendering |
| GEMMA | 0.98.5 | LMM GWAS, kinship matrices, GREML h² (microbiability) |
| PLINK | v1.90b6.21 | Variant QC, LD pruning, PCA, permutation (maxT), allele dose recoding |
| bcftools / samtools | 1.24 / 1.24 | VCF normalization, +setGT / +fill-tags, region extraction |
| R: phyloseq | 1.46.0 | Microbiome data structure, diversity indices |
| R: vegan | 2.7.1 | Richness estimation, Bray–Curtis, PCoA, adonis2 |
| R: randomForest | 4.7.1.2 | Regression/classification of MFD |
| R: lavaan | 0.7.2 | Structural equation modelling (bootstrap 1,000) |
| R: ggplot2 / patchwork | 3.5.2 / 1.3.2 | Graphics and multi-panel assembly |
| R: ggpubr / agricolae / pROC / pdp | 0.6.1 / 1.3.7 / 1.19.0.1 / 0.8.2 | Significance brackets, post-hoc letters, ROC, partial dependence |
| SILVA | Release 138.1 | 16S taxonomy annotation (RDP Classifier 2.2, confidence 0.7) |
| Trimmomatic / FLASH, QIIME 2 / DADA2, PICRUSt2 | performed by the sequencing provider (PE300 workflow) | Read trimming/merging, ASV inference, KEGG functional prediction |

---

## 2. Host genotyping and GWAS

1. **Variant QC and phenotype mapping.** VCF genotypes with depth < 2 or GQ < 10 set to missing; HWE, MAF, allele number, missingness and heterozygosity recalculated (bcftools +setGT / +fill-tags). Loci with > 0.30 missing rate, MAF < 0.05, mean depth outside 4–1000× or significant excess heterozygosity excluded; scaffold contigs dropped, chromosome 1–38 + Z retained; SNP IDs standardized as `chr:pos:ref:alt`. The MFD phenotype was joined to genotypes through the sample mapping table and written into the .fam pedigree column; individuals with both phenotype and genotypes (n = 289) formed the GWAS set.
2. **Effective number of independent tests (Meff).** LD pruning (`--indep-pairwise 50 5 0.2`) yielded 610,751 independent SNPs; genome-wide and suggestive thresholds were set at 0.05/Meff = 8.19 × 10⁻⁸ and 1/Meff = 1.64 × 10⁻⁶. A Li–Ji eigenvalue-based Me calculation was implemented as a cross-check.
3. **Population structure.** Principal components from the LD-pruned genotype matrix (PLINK `--pca`); PC1 (26.1%) and PC2 (12.2%) showed no discrete subgroups and were used as covariates.
4. **GWAS.** Centered GRM built with GEMMA `-gk 1`; genome-wide scan with `-lmm 4` (PC1, PC2 as covariates; 11,235,592 SNPs tested). BH-FDR and Meff-based thresholds applied in R; Manhattan and QQ plots generated from the scan output.
5. **Permutation-based empirical threshold.** PLINK `--assoc` with maxT permutation (1,000 resamples, fixed seed); the 95th percentile of per-permutation genome-wide minimum P values provided the empirical significance threshold.
6. **Regional support.** Pairwise r² among SNPs within ±500 kb of the lead SNP (PLINK `--r2`) rendered as an LD heatmap with gene models and a two-track LocusZoom-style plot; segment-level coverage tracks (CDS/exon/mRNA) extracted from the sorted, tabix-indexed GFF3 annotation.
7. **Lead SNP allele effect.** Additive dosage (0/1/2) of the top SNPs recoded with PLINK `--recode A`, merged with MFD; genotype-class means compared by ANOVA / linear trend and displayed as boxplots.
8. **Suggestive SNP annotation.** Genes within ±5 kb of each suggestive SNP retrieved with tabix and classified as within / upstream / downstream with distance.

## 3. Microbiome (16S) processing and analyses

9. **Microbiome data structure.** ASV table and taxonomy (provider-delivered, QIIME 2/DADA2 inferred; 47,899 ASVs before filtering) imported into phyloseq with a sample metadata table (segment, chicken ID, mapped MFD). ASVs prevalent in < 5% of samples removed (1,042 ASVs retained); relative abundance and CLR transformations derived from counts.
10. **Diversity and spatial structure.** Observed richness and Shannon index compared among segments (Kruskal–Wallis with agricolae post-hoc lettering); Bray–Curtis PCoA with adonis2 PERMANOVA (999 permutations); top-30 annotated ASVs shown as a Z-scored abundance heatmap.
11. **MFD extreme-group contrast.** Birds in the bottom/top 20% of MFD defined as low/high classes; alpha diversity compared per segment (Wilcoxon rank-sum with BH correction over 5 segments × 2 indices), beta diversity by PERMANOVA, and each filtered ASV correlated with MFD (Spearman).
12. **Random forest prediction of MFD.** CLR abundances of all five segments flattened per bird (individuals beyond ±3 SD removed); supervised pre-filtering by Spearman P < 0.05; regression (1,000 trees, node size 5, seed 123, OOB R²/RMSE) with variable importance and partial dependence; classification of high/low MFD quartiles with multi-segment vs cecum-only ROC comparison (pROC).
13. **Microbiability.** Genus × segment traits (prevalence ≥ 30% within segment, mean relative abundance ≥ 0.05%; CLR-transformed, 422 traits) fitted as responses in GEMMA zero models with the GRM to estimate h² (GREML); results summarized as a heatmap.
14. **Microbial GWAS (mGWAS).** Six RF-selected key taxa (ASV-level CLR) scanned against all genotyped SNPs with GEMMA `-lmm 4`; the Meff for the sparse 100k-SNP subset obtained by LD pruning; Bonferroni/BH thresholds applied; bipartite networks linking each taxon's suggestive SNPs (±10 kb) to nearest host genes; lead-SNP dosage effects on taxon CLR abundance displayed as multi-panel boxplots.

## 4. Integration analyses

15. **Mediation.** SNP–taxon–MFD triads tested by the product-of-coefficients framework with custom R code: a path (taxon ~ SNP dosage + PC1 + PC2), b/c′ paths (MFD ~ taxon + SNP dosage), total effect c; indirect effect ab bootstrapped over 2,000 individual resamples for 95% CI and mediation proportion; intersecting GWAS∩mGWAS SNPs verified in a classical mediation model.
16. **Structural equation modelling.** Standardized paths fitted with lavaan (bootstrap 1,000): a latent SCFA-producing module (first PC of four SCFA genera) and single-taxon mediators on the MFD lead-SNP → mediator → MFD chain, with indirect-effect CIs.
17. **Functional association.** Provider-predicted KEGG level-3 pathway abundances (PICRUSt2; 273 pathways) CLR-transformed and regressed on MFD within each segment with PC1/PC2 covariates (1,365 tests), BH-corrected within segment; SCFA/lipid-metabolism candidates visualized across segments together with genus–pathway coupling.

## 5. Figure and table production

18. **Figures.** Multi-panel figures assembled with patchwork/magick: Fig. 1 (alpha/beta diversity, heatmaps, extreme groups), Fig. 2 (random forest), Fig. 3 (GWAS: PCA, Manhattan, LD, allele effect), Fig. 4 (microbiability, mGWAS, networks), Fig. 5 (mediation/SEM), Fig. 6 (functional prediction); supplementary figures S1–S11 (composition, curves, Venn, allele effect, sensitivity). Main figures exported as 300-dpi PDF and 200-dpi PNG; all supplementary tables written as delimited text.

## Reproducibility notes

- Random seeds fixed (e.g., 123 for random forest; fixed seed for PLINK permutations; bootstrap resampling with fixed seeds).
- All P-values two-sided; multiple-testing corrections (BH-FDR, Meff/Bonferroni) applied within the stated test families.
- Segments abbreviated D/K/H/M/R for duodenum/jejunum/ileum/cecum/rectum in raw sample codes.
- Provider-dependent steps (WGS read alignment and variant calling, 16S read processing/DADA2, PICRUSt2) were performed by the sequencing provider; their software versions and the reference genome assembly are documented in the vendor QC report.
