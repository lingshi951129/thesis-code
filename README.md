# thesis-code
This is the repository includes all homemade scripts used in thesis.
This is a short-term public repository which will be turned into private after thesis examination.

For TCGA-LAML analyses, the codes are:
1. assign_gene_Knee_as_ref.pl to assign genes into their domains.
2. innermost.pl to find the inner most TAD for genes that can assign to  multiple domains.
3. determ_domain_condi_Knee_as_ref.pl to determine the same domain and cross boundary gene pairs.
4. correlation_cound_TPM_flt.py to calculate gene correlations for same domain and cross boundary gene pairs for dysreguated boundary identification with the filtering of low expressed genes (more than half of the gene expression is zero).
5. correlation_cound_TPM_wo_flt.py to calculate gene correlations for same domain and cross boundary gene pairs for 500kb gene distance correlation changes analysis without the filtering of low expressed genes (more than half of the gene expression is zero).
6. delta_correlation_wo_flt.pl calculate the delta correlation for 500kb gene distance correlation changes
7. plot_line.R to plot the line plot for 500kb gene distance correlation changes
8. volcano_plot.R to plot the volcano plot 500kb gene distance correlation changes
9. compare_correlation_flt.pl to detect the same domain and cross boundary correlation changes for dysreguated boundary identification analysis
10. dysregulated_TADs_flt.pl to identify dysreguated boundaries.

RNA-Seq analyses code:
1. raw_count_to_FPKM.R to calculate FPKM for htseq counts files
2. FPKM_CVT_TPM.pl to convert FPKM to TPM
3. DEG.R to do the differential gene expression by edgeR
4. boxplot.R to plot the box plot for TCGA-LAML TPM values for DNMT3A mutant and wild type cases.

TAD and loop track making and comparison code:
1. TAD_track.pl for making the UCSC genome browser tracks by Juicer Arrowhead output files
2. loop_track.pl for making the UCSC genome browser tracks by Juicer HiCCUPS output files
3. TAD_compare.pl to compare TADs between two samples
4. loop_compare.pl to compare loops between two samples
