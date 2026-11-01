# =============================================================================
# 04b_DAA_16S_merged_dotplot.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Part A  Build ONE long data frame (one row per OTU x sample) with
#             - read counts and relative abundance
#             - sample metadata (treatment, time, rep, ...)
#             - taxonomy of each OTU (Kingdom ... Species)
#             - differential abundance (DAA) results of DESeq2 AND ANCOM-BC2
#   Part B  Draw a family-level dot plot of the OTUs selected by both methods
#           (family on the y axis, log2 fold change on the x axis, one panel
#           per amendment; the two methods of the same OTU are joined by a line).
#
# BEFORE YOU RUN IT
#   04_DAA_16S.Rmd must have been run once, so that it has SAVED the results of
#   both methods in 03_Analysis_files/ (one file per amendment):
#     deseq_{chicken,cow,compost,pellet}_bact.RDa     (DESeq2)
#     ancombc_{chicken,cow,compost,pellet}_bact.RDa   (ANCOM-BC2)
#   The script stops with a clear message if one of them is missing.
#
# HOW IT RELATES TO 04_DAA_16S.Rmd
#   The Rmd is NOT modified by this script. The script "sources" it: it reads the
#   Rmd and runs only the chunks that load and filter the data and compute the
#   selection thresholds. Chunks that fit models (DESeq2, ANCOM-BC2), draw plots
#   or write files are skipped.
#
# OUTPUT
#   03_Analysis_files/DAA_16S_merged_long.rds / .csv         merged table (Part A)
#   04_Figures/DAA_16S_family_dotplot.png                    figure (Part B)
#   03_Analysis_files/DAA_16S_family_dotplot.rds             the ggplot object
#   03_Analysis_files/DAA_16S_family_dotplot_data.csv        the plotted points
#
# KEY DESIGN CHOICES
#   * DAA results exist per OTU x contrast (amendment vs control), not per
#     sample. Each sample gets the results of ITS OWN treatment's contrast.
#     Control samples have NA in the DAA columns (they are the reference).
#   * Log fold changes are harmonised so that POSITIVE = enriched in the
#     amendment compared with the control, for all four contrasts and both
#     methods. The values exactly as saved by the Rmd are kept in *_lfc_original.
#   * OTUs that could not be tested (NA p-values) are KEPT, with NA results.
#
# HOW TO READ THE FIGURE (Part B)
#   * One panel per amendment. y axis: families with at least one selected OTU,
#     grouped by phylum and ordered by mean log2FC, with the number of selected
#     OTUs in parentheses, e.g. "Bacillaceae (5)".
#   * x axis: log2 fold change (amendment vs control); positive = enriched in
#     the amendment.
#   * One point per OTU and method; only OTUs that PASSED the selection
#     thresholds (adjusted p < 0.05 and |LFC| above the 97.5th percentile, as in
#     the volcano plots) are drawn. By default only OTUs selected by BOTH methods
#     in that contrast are kept (`only_both_methods`).
#   * Shape = method (circle DESeq2, triangle ANCOM-BC2); a line joins the two
#     points of the same OTU (short, flat line = the methods agree).
#   * Size = -log10(adjusted p). Colour = phylum.
# =============================================================================


# ---- 0. Working directory ---------------------------------------------------
# The Rmd uses paths relative to the project root, so make sure we are there.
if (basename(getwd()) == "02_Scripts") setwd("..")
stopifnot("Run this script from the project root (Banana_clean/)" =
            file.exists("02_Scripts/04_DAA_16S.Rmd"))


# ---- 1. Helper: run selected chunks of an Rmd --------------------------------
# Reads an .Rmd line by line, collects the code of each ```{r} chunk and
# evaluates it in `envir` (the global environment by default), in the original
# order. A chunk is SKIPPED if its code matches `skip_pattern` (a regex).
# Individual code lines matching `drop_pattern` (a regex) are removed from the
# chunks that do run (used below to avoid loading a package we do not need).
source_rmd_chunks <- function(rmd_path, skip_pattern, drop_pattern = NULL,
                              envir = globalenv()) {
  lines       <- readLines(rmd_path, warn = FALSE)
  fence_open  <- grepl("^```\\{r", lines)      # a line like  ```{r label}
  fence_close <- grepl("^```\\s*$", lines)     # a line with just  ```
  in_chunk <- FALSE
  code     <- character()
  header   <- ""

  for (ln in seq_along(lines)) {
    if (!in_chunk && fence_open[ln]) {         # chunk starts
      in_chunk <- TRUE
      code     <- character()
      header   <- lines[ln]
      next
    }
    if (in_chunk && fence_close[ln]) {         # chunk ends -> decide: run or skip
      in_chunk <- FALSE
      label <- sub("^```\\{r[ ,]*([^,}]*).*$", "\\1", header)  # chunk name, if any
      if (!nzchar(label)) label <- paste("unnamed chunk ending at line", ln)
      if (grepl(skip_pattern, paste(code, collapse = "\n"))) {
        message("  [skip] ", label)
      } else {
        message("  [run]  ", label)
        if (!is.null(drop_pattern)) code <- code[!grepl(drop_pattern, code)]
        eval(parse(text = code), envir = envir)
      }
      next
    }
    if (in_chunk) code <- c(code, lines[ln])   # inside a chunk: collect the code
  }
  invisible(NULL)
}


# ---- 2. Settings you may want to change (dot plot, Part B) --------------------------------------
# TRUE  = keep only OTUs selected by BOTH methods in that contrast (each OTU then
#         has two points joined by a line).
# FALSE = keep OTUs selected by EITHER method (OTUs selected by one method have a
#         single point and no line).
only_both_methods <- TRUE

# Fewest distinct selected OTUs a family needs (summed over all panels) to be
# shown. 1 = show every family with at least one selected OTU.
min_selected_otus <- 1

# ANCOM-BC2 log fold changes are on the NATURAL log scale (the official ANCOM-BC2
# tutorial describes its plots as "log fold-change (in natural log)"). DESeq2 is
# on log2, so we divide the ANCOM-BC2 values by log(2) to put both on log2.
# Set to FALSE only if you ever feed in values that are already log2.
# NB: this only affects where ANCOM-BC2 points are DRAWN. Which OTUs are
# "selected" was decided earlier, on the original scale, and is not changed.
ancombc_is_natural_log <- TRUE

# Vertical layout inside one family row (row height = 1)
method_dodge <- 0   # DESeq2 is drawn +this above, ANCOM-BC2 -this below the centre
slot_span    <- 0.25   # OTUs of the same family spread over +/- this around the centre

# Order and labels of the panels
contrast_levels <- c("chicken", "cow", "compost", "pellet")
contrast_labels <- c(chicken = "Chicken manure", cow = "Cow manure",
                     compost = "Compost",       pellet = "Pelletised compost")

# Phylum colours: the same palette used for the taxonomy barplot (plotColorBarplot_16S).
# "Chloroflexota" (the name in the dada2 taxonomy) shares the colour of
# "Chloroflexi", and "Unclassified" (OTUs with no phylum) uses the "Unassigned" grey.
phylum_colours <- c(
  "Cyanobacteria"     = "#8DD3C7",
  "Actinomycetota"    = "#B2D8EE",
  "Bacteroidota"      = "#BEBADA",
  "Acidobacteriota"   = "#FB8072",
  "Bacillota"         = "#80B1D3",
  "Gemmatimonadota"   = "#FDB462",
  "Planctomycetota"   = "#B3DE69",
  "Pseudomonadota"    = "#FCCDE5",
  "Myxococcota"       = "#256D3D",
  "Chloroflexi"       = "#A34022",
  "Chloroflexota"     = "#A34022",
  "Verrucomicrobiota" = "#7C4D79",
  "Unclassified"      = "gray30"
)
# Phyla in the figure that are NOT in the palette above get one of these, in
# order (the first is the "Less than 2%" colour of the barplot). A message
# lists which phyla used a fallback.
fallback_colours <- c("#F9E77F", "#C9B79C", "#7F7F7F")

out_png  <- "04_Figures/DAA_16S_family_dotplot.png"
out_plot <- "03_Analysis_files/DAA_16S_family_dotplot.rds"
out_data <- "03_Analysis_files/DAA_16S_family_dotplot_data.csv"


# ---- 3. Load the saved DAA results -------------------------------------------
# 04_DAA_16S.Rmd saves the results of each method, one file per amendment:
#   DESeq2    deseq_*_bact.RDa    data.frame: OTU, baseMean, log2FoldChange, lfcSE,
#                                 pvalue, padj, treatment
#   ANCOM-BC2 ancombc_*_bact.RDa  data.frame: taxon, lfc, se, p_val, p_adj, treatment
# The objects are given the names the Rmd uses (res_chicken, ... and chicken, ...)
# because the "deseq bind" and "ancombc bind" chunks that we run below need them.
amendments  <- c("chicken", "cow", "compost", "pellet")
saved_files <- c(sprintf("03_Analysis_files/deseq_%s_bact.RDa",   amendments),
                 sprintf("03_Analysis_files/ancombc_%s_bact.RDa", amendments))
missing_files <- saved_files[!file.exists(saved_files)]
if (length(missing_files) > 0)
  stop("Saved results not found:\n  ", paste(missing_files, collapse = "\n  "),
       "\nRun the DESeq2 / ANCOM-BC2 chunks of 02_Scripts/04_DAA_16S.Rmd first.")

res_chicken <- readRDS("03_Analysis_files/deseq_chicken_bact.RDa")
res_cow     <- readRDS("03_Analysis_files/deseq_cow_bact.RDa")
res_compost <- readRDS("03_Analysis_files/deseq_compost_bact.RDa")
res_pellet  <- readRDS("03_Analysis_files/deseq_pellet_bact.RDa")

chicken <- readRDS("03_Analysis_files/ancombc_chicken_bact.RDa")
cow     <- readRDS("03_Analysis_files/ancombc_cow_bact.RDa")
compost <- readRDS("03_Analysis_files/ancombc_compost_bact.RDa")
pellet  <- readRDS("03_Analysis_files/ancombc_pellet_bact.RDa")


# ---- 4. Source the original Rmd (selected chunks only) -----------------------
# Skipped = any chunk that contains one of:
#   DESeq( / DESeqDataSetFromMatrix   -> fitting the DESeq2 models (results are saved)
#   ancombc2(                         -> fitting ANCOM-BC2 (slow; results are saved)
#   ggplot( / ggsave                  -> volcano plots and writing figures
#   saveRDS                           -> overwriting saved results / figures
#   readRDS("~                        -> figure objects read from an absolute path
# What still runs: libraries, colours, data loading, taxa filtering, set.seed,
# the sample subsets, and the two "bind" chunks (which define the fold-change
# thresholds `threshold_deseq` and `threshold_ancombc`).
skip_pattern <- paste(c("DESeq\\(", "DESeqDataSetFromMatrix", "ancombc2\\(",
                        "ggplot\\(", "ggsave", "saveRDS", "readRDS\\(\"~"),
                      collapse = "|")

# Package lines we do not need (nothing is fitted here) are dropped from the
# chunks that run: this is faster and avoids version conflicts (e.g. ANCOMBC/CVXR)
# and function-name masking by the Bioconductor packages.
drop_pattern <- "^\\s*library\\((ANCOMBC|DESeq2|apeglm)\\)"

message("Sourcing 02_Scripts/04_DAA_16S.Rmd ...")
source_rmd_chunks("02_Scripts/04_DAA_16S.Rmd", skip_pattern, drop_pattern)

# Objects now available in the session (created by the Rmd):
#   zotu97                      filtered phyloseq object (counts, taxonomy, sample data)
#   treat.order, treat.colors   treatment order / colours
#   threshold_deseq / threshold_ancombc   |LFC| cut-offs (97.5th percentile)
#   number_selected_deseq / number_selected_ancombc   OTUs above the thresholds
library(tidyverse)   # (already loaded by the Rmd; repeated for clarity)
library(phyloseq)


# ---- 5. DESeq2 results -> one tidy table (OTU x contrast) --------------------
# Sign convention. We want POSITIVE = enriched in the amendment (amendment vs
# control) for every contrast.
#  * Files saved by the CURRENT 04_DAA_16S.Rmd (control = reference level) already
#    have this direction and carry the tag attr(, "lfc_direction") =
#    "amendment_vs_control": no flip (factor +1).
#  * Files saved by the OLD version of the Rmd (no tag) used the alphabetical
#    reference level: for chicken and compost the coefficient was "control vs
#    amendment", i.e. REVERSED (checked against raw counts: those OTUs were
#    more abundant in the amendment). For those we multiply by -1 (legacy flip).
# NB: using `get_flip()` also protects against flipping twice.
get_flip <- function(df, legacy_flip) {
  if (identical(attr(df, "lfc_direction"), "amendment_vs_control")) 1 else legacy_flip
}
deseq_flip <- c(chicken = get_flip(res_chicken, -1), cow     = get_flip(res_cow,     1),
                compost = get_flip(res_compost, -1), pellet  = get_flip(res_pellet,  1))

deseq_tidy <- bind_rows(res_chicken, res_cow, res_compost, res_pellet) %>%
  # NB: unlike the Rmd, we do NOT na.omit() here, so untestable OTUs are kept
  transmute(
    OTU,
    contrast             = treatment,                           # amendment name
    deseq_baseMean       = baseMean,                            # mean normalised count
    deseq_lfc_original   = log2FoldChange,                      # as in the Rmd
    deseq_lfc            = log2FoldChange * unname(deseq_flip[treatment]),  # + = up in amendment
    deseq_se             = lfcSE,
    deseq_pvalue         = pvalue,
    deseq_padj           = padj,
    # same "selected" rule as the Rmd volcano plot: padj < 0.05 and |LFC| above
    # the 97.5th-percentile threshold (the sign does not matter for |LFC|)
    deseq_selected       = ifelse(padj < 0.05 & abs(log2FoldChange) > threshold_deseq,
                                  treatment, "no")
  )


# ---- 6. ANCOM-BC2 results -> one tidy table (OTU x contrast) -----------------
# Same sign logic as for DESeq2 (see above): files saved by the current Rmd are
# tagged "amendment_vs_control" (no flip); older chicken and compost files held
# `lfc_treatmentcontrol` = control vs amendment and get the legacy flip (-1).
ancombc_flip <- c(chicken = get_flip(chicken, -1), cow     = get_flip(cow,     1),
                  compost = get_flip(compost, -1), pellet  = get_flip(pellet,  1))

# report the sign factors used (-1 = untagged file from the OLD Rmd, flipped here)
message("DESeq2 sign factors:    ", paste(names(deseq_flip),   deseq_flip,   collapse = ", "))
message("ANCOM-BC2 sign factors: ", paste(names(ancombc_flip), ancombc_flip, collapse = ", "))

ancombc_tidy <- bind_rows(chicken, cow, compost, pellet) %>%
  transmute(
    OTU                  = taxon,
    contrast             = treatment,
    ancombc_lfc_original = lfc,                                 # as saved by the Rmd
    ancombc_lfc          = lfc * unname(ancombc_flip[treatment]),  # + = up in amendment
    ancombc_se           = se,
    ancombc_pvalue       = p_val,
    ancombc_padj         = p_adj,
    ancombc_selected     = ifelse(p_adj < 0.05 & abs(lfc) > threshold_ancombc,
                                  treatment, "no")
  )

# Sanity check: our "selected" counts must match the counts printed in the Rmd
# plot subtitles (number_selected_deseq / number_selected_ancombc).
if (sum(deseq_tidy$deseq_selected != "no", na.rm = TRUE) != number_selected_deseq)
  warning("DESeq2 'selected' count differs from the Rmd")
if (sum(ancombc_tidy$ancombc_selected != "no", na.rm = TRUE) != number_selected_ancombc)
  warning("ANCOM-BC2 'selected' count differs from the Rmd")


# ---- 7. Combine both methods (OTU x contrast) --------------------------------
# full_join: keeps an OTU/contrast even if only one method tested it.
daa <- full_join(deseq_tidy, ancombc_tidy, by = c("OTU", "contrast")) %>%
  mutate(
    # which method(s) flagged the OTU as differentially abundant?
    selected_by = case_when(
      deseq_selected %in% treat.order[-1] & ancombc_selected %in% treat.order[-1] ~ "both",
      deseq_selected %in% treat.order[-1]                                          ~ "deseq2_only",
      ancombc_selected %in% treat.order[-1]                                        ~ "ancombc2_only",
      TRUE                                                                         ~ "none"
    )
  )
stopifnot("Duplicated OTU x contrast rows in DAA table" =
            !anyDuplicated(daa[, c("OTU", "contrast")]))


# ---- 8. Counts, sample data and taxonomy from the phyloseq object ------------
# 7a. Counts: OTU x sample matrix -> long format (one row per OTU x sample).
otu_mat <- as(otu_table(zotu97), "matrix")
if (!taxa_are_rows(zotu97)) otu_mat <- t(otu_mat)    # make sure rows = OTUs

counts_long <- otu_mat %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  # relative abundance within the retained (filtered) OTUs of each sample
  mutate(rel_abundance = count / sum(count)) %>%
  ungroup()

# 7b. Sample metadata (treatment, time, rep, sample_ID, treat.pot).
sample_df <- data.frame(sample_data(zotu97)) %>%
  rownames_to_column("sample")

# 7c. Taxonomy (dada2 assignment loaded in the Rmd) + the lowest assigned rank.
tax_df <- as(tax_table(zotu97), "matrix") %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
  mutate(lowest_rank = coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom))


# ---- 9. Assemble the final table ---------------------------------------------
# The count table contains samples without metadata (negative / mock controls
# and samples outside this experiment). inner_join() keeps only samples that have
# metadata, i.e. the experimental samples used in the DAA.
merged_16S <- counts_long %>%
  inner_join(sample_df, by = "sample") %>%
  dplyr::filter(!is.na(treatment)) %>%
  left_join(tax_df, by = "OTU") %>%
  # Each sample gets the DAA results of its own contrast. Controls are the
  # reference, so they have no contrast (NA) and therefore NA in the DAA columns.
  mutate(contrast = if_else(treatment == "control", NA_character_, as.character(treatment))) %>%
  left_join(daa, by = c("OTU", "contrast"), na_matches = "never") %>%
  # tidy column order: identifiers & counts | sample data | taxonomy | DAA
  relocate(OTU, sample, count, rel_abundance,
           any_of(c("treatment", "time", "rep", "sample_ID", "treat.pot")),
           any_of(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                    "Species", "lowest_rank")),
           contrast)


# ---- 10. Checks ---------------------------------------------------------------
# (a) one row per OTU x sample, nothing lost or duplicated by the joins
n_otu    <- n_distinct(counts_long$OTU)
n_sample <- n_distinct(merged_16S$sample)
stopifnot("Row count != OTUs x samples" = nrow(merged_16S) == n_otu * n_sample,
          "Duplicated OTU x sample rows" = !anyDuplicated(merged_16S[, c("OTU", "sample")]))
message(sprintf("Final table: %d rows = %d OTUs x %d samples, %d columns",
                nrow(merged_16S), n_otu, n_sample, ncol(merged_16S)))

# (b) independent check of the fold-change DIRECTION. For each OTU and contrast,
# compute the empirical log2 ratio of mean relative abundance (amendment /
# control) straight from the counts, and correlate it with each method's
# harmonised LFC. Correlations should be POSITIVE for every contrast.
mean_ra <- merged_16S %>%
  group_by(OTU, treatment) %>%
  summarise(m = mean(rel_abundance), .groups = "drop") %>%
  pivot_wider(names_from = treatment, values_from = m)

direction_check <- mean_ra %>%
  pivot_longer(all_of(setdiff(treat.order, "control")),
               names_to = "contrast", values_to = "m_amendment") %>%
  mutate(empirical_log2fc = log2((m_amendment + 1e-6) / (control + 1e-6))) %>%
  inner_join(daa, by = c("OTU", "contrast")) %>%
  group_by(contrast) %>%
  summarise(spearman_deseq   = cor(deseq_lfc,   empirical_log2fc, method = "spearman", use = "complete.obs"),
            spearman_ancombc = cor(ancombc_lfc, empirical_log2fc, method = "spearman", use = "complete.obs"),
            .groups = "drop")
print(direction_check)
if (any(direction_check[, -1] < 0, na.rm = TRUE))
  warning("Negative correlation with the empirical fold change: check LFC signs!")


# ---- 11. Save ----------------------------------------------------------------
saveRDS(merged_16S, "03_Analysis_files/DAA_16S_merged_long.rds")
write_csv(merged_16S, "03_Analysis_files/DAA_16S_merged_long.csv")
message("Saved to 03_Analysis_files/DAA_16S_merged_long.{rds,csv}")

# =============================================================================
# PART B - Family-level dot plot of the OTUs selected by the DAA methods
# =============================================================================


# ---- 12. One row per OTU x contrast -------------------------------------------
# The merged table has one row per OTU x SAMPLE; the DAA columns are repeated on
# every sample of a contrast, so we collapse to unique OTU x contrast rows.
# Control samples have no contrast (NA) and are dropped here.
otu_contrast <- merged_16S %>%
  dplyr::filter(!is.na(contrast)) %>%
  distinct(OTU, contrast, Phylum, Family, lowest_rank,
           deseq_lfc,   deseq_padj,   deseq_selected,
           ancombc_lfc, ancombc_padj, ancombc_selected)
stopifnot("OTU x contrast is not unique" = !anyDuplicated(otu_contrast[, c("OTU", "contrast")]))


# ---- 13. Long format: one row per OTU x contrast x METHOD ---------------------
# Keep a point only if THAT method selected the OTU in THAT contrast.
lfc_scale_ancombc <- if (ancombc_is_natural_log) log(2) else 1

points_all <- bind_rows(
  otu_contrast %>% transmute(OTU, contrast, Phylum, Family, lowest_rank,
                    method = "DESeq2",
                    lfc    = deseq_lfc,
                    padj   = deseq_padj,
                    selected = deseq_selected != "no"),      # NA if not tested
  otu_contrast %>% transmute(OTU, contrast, Phylum, Family, lowest_rank,
                    method = "ANCOM-BC2",
                    lfc    = ancombc_lfc / lfc_scale_ancombc,   # -> log2 scale
                    padj   = ancombc_padj,
                    selected = ancombc_selected != "no")
) %>%
  dplyr::filter(selected %in% TRUE) %>%          # drops FALSE and NA
  mutate(
    method   = factor(method, levels = c("DESeq2", "ANCOM-BC2")),
    contrast = factor(contrast, levels = contrast_levels),
    Phylum   = coalesce(Phylum, "Unclassified"),
    # Row label = family. OTUs without a family get "<lowest assigned rank>
    # (unclass.)", e.g. "Rhizobiales (unclass.)", instead of being lumped together.
    taxon_label = coalesce(Family, paste0(lowest_rank, " (unclass.)")),
    neg_log_padj = -log10(pmax(padj, 1e-300))     # guard against padj == 0
  )

# optionally keep only OTUs selected by both methods (same OTU, same contrast)
if (only_both_methods) {
  points_all <- points_all %>%
    group_by(OTU, contrast) %>%
    dplyr::filter(n_distinct(method) == 2) %>%
    ungroup()
}

# keep only families with enough selected OTUs
keep_families <- points_all %>%
  group_by(taxon_label) %>%
  summarise(n_otu = n_distinct(OTU), .groups = "drop") %>%
  dplyr::filter(n_otu >= min_selected_otus) %>%
  pull(taxon_label)
points_all <- dplyr::filter(points_all, taxon_label %in% keep_families)


# ---- 14. Row order: grouped by phylum, then by mean log2FC --------------------
# Phyla ordered by number of selected OTUs (most frequent first); the rows and
# the legend follow this order.
phylum_levels <- points_all %>%
  distinct(OTU, Phylum) %>%
  dplyr::count(Phylum, sort = TRUE) %>%
  pull(Phylum)

family_order <- points_all %>%
  group_by(taxon_label) %>%
  summarise(Phylum   = dplyr::first(Phylum),
            mean_lfc = mean(lfc),
            n_otu    = n_distinct(OTU),     # distinct selected OTUs of the family
            .groups  = "drop") %>%
  mutate(phylum_rank = match(Phylum, phylum_levels),
         axis_label  = paste0(taxon_label, " (", n_otu, ")")) %>%   # e.g. "Bacillaceae (5)"
  arrange(phylum_rank, mean_lfc) %>%
  # y position of the family row: first row on top (largest y value)
  mutate(y_base = rev(seq_len(n())))


# ---- 15. Colour groups --------------------------------------------------------
# one colour per phylum: from the barplot palette where defined, else a fallback
missing_phyla <- setdiff(phylum_levels, names(phylum_colours))
stopifnot("More unlisted phyla than `fallback_colours`" =
            length(missing_phyla) <= length(fallback_colours))
if (length(missing_phyla) > 0)
  message("No palette colour for: ", paste(missing_phyla, collapse = ", "),
          " -> using fallback colours")
phylum_palette <- c(phylum_colours[intersect(phylum_levels, names(phylum_colours))],
                    setNames(fallback_colours[seq_along(missing_phyla)], missing_phyla))
phylum_palette <- phylum_palette[phylum_levels]     # legend order = frequency order

points_all <- points_all %>%
  mutate(phylum_group = factor(Phylum, levels = phylum_levels))


# ---- 16. Vertical position of every point -------------------------------------
# (a) OTU slot: OTUs of the same family in the same panel are spread evenly over
#     +/- slot_span, ordered by their mean log2FC (lowest at the bottom).
otu_slots <- points_all %>%
  group_by(taxon_label, contrast, OTU) %>%
  summarise(otu_lfc = mean(lfc), .groups = "drop") %>%
  group_by(taxon_label, contrast) %>%
  arrange(otu_lfc, .by_group = TRUE) %>%
  mutate(n_otu      = n(),
         slot       = row_number(),
         otu_offset = if_else(n_otu == 1, 0,
                              -slot_span + (slot - 1) * 2 * slot_span / (n_otu - 1))) %>%
  ungroup() %>%
  dplyr::select(taxon_label, contrast, OTU, otu_offset)

# (b) method dodge: DESeq2 slightly above, ANCOM-BC2 slightly below the OTU slot.
plot_df <- points_all %>%
  left_join(family_order[, c("taxon_label", "y_base")], by = "taxon_label") %>%
  left_join(otu_slots, by = c("taxon_label", "contrast", "OTU")) %>%
  mutate(y = y_base + otu_offset + if_else(method == "DESeq2", method_dodge, -method_dodge))


# ---- 17. Lines joining the two methods (the same OTU, selected by both) -------
link_df <- plot_df %>%
  mutate(m = if_else(method == "DESeq2", "deseq", "ancombc")) %>%
  dplyr::select(OTU, contrast, m, lfc, y) %>%
  pivot_wider(names_from = m, values_from = c(lfc, y)) %>%
  drop_na()          # keeps only OTUs that have a point for both methods


# ---- 18. Plot -----------------------------------------------------------------
dotplot <- ggplot(plot_df, aes(x = lfc, y = y)) +
  # zero line: left of it = depleted, right of it = enriched in the amendment
  geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.5) +
  # line between the two methods for OTUs selected by both (drawn first = behind)
  geom_segment(data = link_df,
               aes(x = lfc_deseq, xend = lfc_ancombc, y = y_deseq, yend = y_ancombc),
               inherit.aes = FALSE, colour = "#52514e", linewidth = 0.4) +
  # points: fill = phylum, shape = method, size = significance.
  # shapes 21/24 are fillable; the thin dark outline keeps the pastel fills visible
  # on a white background and separates overlapping points
  geom_point(aes(size = neg_log_padj, shape = method, fill = phylum_group),
             colour = "#52514e", stroke = 0.4, alpha = 0.95) +
  facet_grid(. ~ contrast, labeller = labeller(contrast = contrast_labels)) +
  scale_shape_manual(values = c("DESeq2" = 21, "ANCOM-BC2" = 24), name = "Method") +
  scale_fill_manual(values = phylum_palette, name = "Phylum", drop = FALSE) +
  scale_size_continuous(range = c(1.8, 6), trans = "sqrt", breaks = c(2, 5, 10, 25, 50, 100), name = expression(-log[10]~"(adj. p)")) +
  scale_y_continuous(breaks = family_order$y_base, labels = family_order$axis_label,
                     expand = expansion(add = 0.7)) +
  # legends: make shapes visible in the colour legend and vice versa
  guides(fill  = guide_legend(order = 1, override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(order = 2, override.aes = list(fill = "#898781", size = 4)),
         size  = guide_legend(order = 3)) +
  labs(x = "log2 fold change (amendment vs control)", y = "Family") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_line(colour = "#e1e0d9", linewidth = 0.3),
    panel.grid.major.y  = element_line(colour = "#c3c2b7", linewidth = 0.3, linetype = "dashed"),
    panel.border        = element_rect(colour = "#c3c2b7", fill = NA),
    axis.text.y         = element_text(face = "italic", colour = "#0b0b0b"),
    axis.text.x         = element_text(colour = "#52514e"),
    strip.text          = element_text(size = 12, face = "bold", colour = "#0b0b0b"),
    plot.caption        = element_text(colour = "#52514e", hjust = 0),
    legend.position     = "right"
  )

dotplot


# ---- 19. Save -----------------------------------------------------------------
n_fam <- nrow(family_order)
ggsave(out_png, dotplot, width = 13, height = max(5, 0.22 * n_fam + 2.5),
       dpi = 300, bg = "white")        # taller figure when there are more families
saveRDS(dotplot, out_plot)
write_csv(dplyr::select(plot_df, OTU, contrast, method, taxon_label, Phylum,
                        lfc, padj, neg_log_padj, y),
          out_data)

message(sprintf("%d families, %d OTU x contrast hits, %d of them selected by both methods",
                n_fam, n_distinct(plot_df$OTU, plot_df$contrast), nrow(link_df)))
