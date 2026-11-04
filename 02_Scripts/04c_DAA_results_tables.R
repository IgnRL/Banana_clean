# =============================================================================
# 04c_DAA_results_tables.R
# Word tables with the differential-abundance results (DESeq2 + ANCOM-BC2)
# -----------------------------------------------------------------------------
# One table per dataset (16S and ITS), exported to Word:
#   04_Tables/DAA_16S_table.docx
#   04_Tables/DAA_ITS_table.docx
#
# Columns: OTU | Phylum | Taxonomy | Treatment | DESeq2 (log2FC, FDR) | ANCOMBC (log2FC, FDR)
#   * Taxonomy = lowest assigned taxon and its rank, e.g. "Mitochondria (Family)".
#     Species names are written in italics ("Aspergillus fumigatus" (Species)).
#   * An OTU that appears for several treatments is shown once: its OTU, Phylum
#     and Taxonomy cells are merged over its rows.
#   * log2FC > 0 = enriched in the amendment compared with the control
#     (the harmonised signs from 04b_*_merged_dotplot.R).
#   * ANCOM-BC2 LFCs are on the natural-log scale in the saved results, so they
#     are divided by log(2) to give log2FC (same as in the dot plots).
#   * FDR = adjusted p-value; values below 0.001 are written "<0.001".
#   * Rows are ordered by OTU number (ZOTU2 < ZOTU10 < ZOTU100), then treatment.
#
# REQUIRES: the merged tables made by 04b_DAA_16S_merged_dotplot.R and
#   04b_DAA_ITS_merged_dotplot.R (03_Analysis_files/DAA_{16S,ITS}_merged_long.rds)
# =============================================================================

library(dplyr)
library(stringr)
library(flextable)
library(officer)

if (basename(getwd()) == "02_Scripts") setwd("..")

# ---- Settings ---------------------------------------------------------------
# TRUE  = only OTUs selected by BOTH methods in that treatment (as in the dot plots)
# FALSE = OTUs selected by EITHER method (the other method's values are still shown)
only_both_methods <- FALSE

# ANCOM-BC2 values are natural log -> divide by log(2). FALSE if already log2.
ancombc_is_natural_log <- TRUE

treat_order <- c("chicken", "cow", "compost", "pellet")
out_dir     <- "04_Tables"
dir.create(out_dir, showWarnings = FALSE)

# ---- Helpers ----------------------------------------------------------------
fmt_lfc <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))
fmt_p   <- function(p) ifelse(is.na(p), "NA",
                              ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
is_sel  <- function(x) !is.na(x) & x != "no"      # "no" / NA = not selected

# ---- Build the table for one dataset ----------------------------------------
build_table <- function(tag) {
  merged <- readRDS(sprintf("03_Analysis_files/DAA_%s_merged_long.rds", tag))

  # the DAA columns are repeated on every sample: collapse to OTU x treatment
  daa <- merged %>%
    dplyr::filter(!is.na(contrast)) %>%
    distinct(OTU, contrast, Phylum, Genus, Species, Family, Order, Class,
             lowest_rank, deseq_lfc, deseq_padj, deseq_selected,
             ancombc_lfc, ancombc_padj, ancombc_selected)
  stopifnot("OTU x contrast is not unique" = !anyDuplicated(daa[, c("OTU", "contrast")]))

  sel_deseq   <- is_sel(daa$deseq_selected)
  sel_ancombc <- is_sel(daa$ancombc_selected)
  daa <- daa[if (only_both_methods) sel_deseq & sel_ancombc else sel_deseq | sel_ancombc, ]

  scale_ancombc <- if (ancombc_is_natural_log) log(2) else 1

  daa %>%
    mutate(
      # lowest assigned taxon and its rank (species shown as "Genus species")
      tax_rank = case_when(!is.na(Species) ~ "Species", !is.na(Genus) ~ "Genus",
                        !is.na(Family)  ~ "Family",  !is.na(Order) ~ "Order",
                        !is.na(Class)   ~ "Class",   TRUE ~ "Phylum"),
      tax_name = case_when(tax_rank == "Species" & !is.na(Genus) ~ paste(Genus, Species),
                        TRUE ~ as.character(lowest_rank)),
      is_species = tax_rank == "Species",
      Phylum     = coalesce(Phylum, "Unclassified"),
      otu_num    = as.integer(str_extract(OTU, "\\d+")),
      treat_i    = match(contrast, treat_order)
    ) %>%
    arrange(otu_num, treat_i) %>%
    transmute(OTU, Phylum,
              tax_name, tax_rank, is_species,                       # helper columns (not shown as such)
              taxonomy    = paste0(tax_name, " (", tax_rank, ")"),
              treatment   = contrast,
              deseq_lfc   = fmt_lfc(deseq_lfc),
              deseq_fdr   = fmt_p(deseq_padj),
              ancombc_lfc = fmt_lfc(ancombc_lfc / scale_ancombc),
              ancombc_fdr = fmt_p(ancombc_padj))
}

# ---- Format and export to Word ----------------------------------------------
export_word <- function(tab, tag) {
  shown <- c("OTU", "Phylum", "taxonomy", "treatment",
             "deseq_lfc", "deseq_fdr", "ancombc_lfc", "ancombc_fdr")

  # last row of each OTU -> thin line between OTUs
  last_of_otu <- which(!duplicated(tab$OTU, fromLast = TRUE))
  last_of_otu <- last_of_otu[last_of_otu < nrow(tab)]

  ft <- flextable(tab, col_keys = shown) %>%
    set_header_labels(OTU = "OTU", Phylum = "Phylum", taxonomy = "Taxonomy",
                      treatment = "Treatment",
                      deseq_lfc = "log2FC", deseq_fdr = "FDR",
                      ancombc_lfc = "log2FC", ancombc_fdr = "FDR") %>%
    # top header row: grouped method names
    add_header_row(values = c("OTU", "Phylum", "Taxonomy", "Treatment", "DESeq2", "ANCOMBC"),
                   colwidths = c(1, 1, 1, 1, 2, 2)) %>%
    merge_v(part = "header", j = 1:4) %>%
    # log2FC with a real subscript 2
    compose(part = "header", i = 2, j = c(5, 7),
            value = as_paragraph("log", as_sub("2"), "FC")) %>%
    # species names in italics (the rank in brackets stays upright)
    compose(i = which(tab$is_species), j = "taxonomy",
            value = as_paragraph(as_i(tax_name), " (", tax_rank, ")")) %>%
    # one cell per OTU: merge OTU, Phylum and Taxonomy over the OTU's rows
    merge_v(j = "OTU", target = c("OTU", "Phylum", "taxonomy")) %>%
    theme_booktabs() %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 9, part = "all") %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    align(j = c(1, 4:8), align = "center", part = "body") %>%
    align(j = 2:3, align = "left", part = "body") %>%
    valign(j = 1:3, valign = "center", part = "body") %>%
    hline(i = last_of_otu, border = fp_border(color = "grey70", width = 0.5), part = "body") %>%
    vline(j = c(4, 6), border = fp_border(color = "grey70", width = 0.5), part = "all") %>%
    set_table_properties(layout = "autofit", width = 1)

  doc <- read_docx() %>%
    body_end_section_landscape() %>%
    body_add_flextable(ft)

  out_file <- file.path(out_dir, sprintf("DAA_%s_table.docx", tag))
  print(doc, target = out_file)
  message(sprintf("%s: %d rows (%d OTUs) -> %s", tag, nrow(tab), n_distinct(tab$OTU), out_file))
  invisible(ft)
}

# ---- Run --------------------------------------------------------------------
tab_16S <- build_table("16S")
tab_ITS <- build_table("ITS")

export_word(tab_16S, "16S")
export_word(tab_ITS, "ITS")

print(head(dplyr::select(tab_16S, OTU, Phylum, taxonomy, treatment,
                         deseq_lfc, deseq_fdr, ancombc_lfc, ancombc_fdr), 10))
