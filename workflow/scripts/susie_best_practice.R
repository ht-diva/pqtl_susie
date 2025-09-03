
# =============================
# SuSiE Fine-mapping Input Loader
# =============================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(data.table)  # For fast I/O
  library(susieR)
})


# ---------- User Inputs ----------
proj_path <- "/scratch/dariush.ghasemi/projects/pqtl_susie/"

#path_sumstat <- glue(proj_path, "test/results/fm/tmp/seq.8221.19_22_24234172_24401503_sumstat.csv")
#path_dosage  <- glue(proj_path, "test/results/fm/tmp/seq.8221.19_22_24234172_24401503_dosage.tsv")
path_sumstat <- snakemake@input[["sumstat"]]
path_dosage  <- snakemake@input[["dosage"]]


# Load parameters for susieR model
label_chr <- snakemake@params[["chrcol"]]
susie_min_abs_cor <- snakemake@params[["min_abs_corr"]]
susie_iter <- snakemake@params[["iter"]]
#susie_iter <- 1000
#susie_min_abs_cor <- 0.0

out_cs_summary <- snakemake@output[["cs_summary"]]
out_cs_list <- snakemake@output[["cs_list"]]
out_cs_rds <- snakemake@output[["cs_rds"]]

# ---------- Helper Functions ----------
err_handling <- function(e) { stop("❌ SuSiE failed: ", e$message) }

check_file <- function(path, min_size = 1e3) {
  if (!file.exists(path)) {
    stop(paste("❌ File does not exist:", path))
  }
  size <- file.size(path)
  if (is.na(size) || size < min_size) {
    stop(paste("❌ File is empty or too small:", path))
  }
  message(paste("✅ File exists and size =", round(size/1024, 2), "KB:", path))
  return(TRUE)
}


flip_alleles <- function(beta) {
  return(-beta)
}

is_strand_ambiguous <- function(a1, a2) {
  ambig <- c("AT", "TA", "GC", "CG")
  return(paste0(a1, a2) %in% ambig)
}


check_file(path_sumstat)
check_file(path_dosage)

# ---------- Load Data ----------
# Use fread with explicit arguments to avoid surprises
sumstat <- tryCatch({
  fread(path_sumstat, header = TRUE, sep = "\t", data.table = FALSE)
}, error = function(e) {
  stop("❌ Failed to read sumstat file: ", e$message)
})

dosage <- tryCatch({
  fread(path_dosage, header = TRUE, sep = " ", data.table = FALSE)
}, error = function(e) {
  stop("❌ Failed to read dosage file: ", e$message)
})

# ---------- Basic QC ----------
# rename colum name
colnames(sumstat)[which(names(sumstat) == label_chr)] <- "CHR"

# Check mandatory columns in summary stats
required_sumstat_cols <- c("SNPID", "CHR", "POS", "EA", "NEA", "BETA", "SE", "MLOG10P")
missing_cols <- setdiff(required_sumstat_cols, colnames(sumstat))

# Remove effective allele from column names
dosage <- dosage %>% dplyr::rename_with(~ str_remove(., "_[ATCG]+$"), matches("_[ATCG]+$"))

# Store variants list for the match; remove FID, IID, ..., PHENOTYPE columns
genotype_variants <- names(dosage)[7:ncol(dosage)]

if (length(missing_cols) > 0) {
  stop("❌ Missing required columns in sumstat file: ", paste(missing_cols, collapse = ", "))
}

message("✅ Summary stats and dosage files loaded successfully.")

# ---------- Variant Matching (to avoid allele mismatch) ----------
common_snps <- intersect(sumstat$SNPID, genotype_variants)

if (length(common_snps) == 0) {
  stop("❌ No overlapping SNPs between sumstat and dosage files.")
}

message("✅ Overlapping SNPs found: ", length(common_snps))

# Optional: subset both datasets to common SNPs
sumstat <- sumstat[sumstat$SNPID %in% common_snps, ]
#X <- dosage[, common_snps]

message("✅ Subsetted to common SNPs. Ready for SuSiE.")

# ---------- Allele Alignment ----------

# Ensure alleles are aligned between summary stats and dosage
# Assume dosage file has A1 (effect) and A2 (other) columns if available
#if (!all(c("EA", "A2") %in% colnames(dosage))) {
#  stop("❌ Dosage file must contain 'A1' and 'A2' columns for allele matching.")
#}

# Merge datasets for alignment
#merged <- merge(sumstat, dosage[, c("SNP", "A1", "A2")], by = "SNP", suffixes = c("_sum", "_dos"))

# Drop strand ambiguous SNPs
#merged <- merged[!is_strand_ambiguous(merged$A1_sum, merged$A2_sum), ]
#message("✅ Removed strand ambiguous SNPs. Remaining SNPs: ", nrow(merged))

# Align effect alleles
#flip_idx <- which(merged$A1_sum != merged$A1_dos & merged$A1_sum == merged$A2_dos)

#if (length(flip_idx) > 0) {
#  merged$BETA[flip_idx] <- flip_alleles(merged$BETA[flip_idx])
#  message("✅ Flipped effect sizes for ", length(flip_idx), " SNPs to match dosage alleles.")
#}


# Subset dosage matrix to aligned SNPs only
#dosage <- dosage[dosage$SNP %in% merged$SNP, ]
#merged <- merged[order(match(merged$SNP, dosage$SNP)), ]
#dosage <- dosage[order(match(dosage$SNP, merged$SNP)), ]

#stopifnot(all(merged$SNP == dosage$SNP)) # safety check

# ---------- Prepare Inputs for SuSiE ----------
betas    <- sumstat$BETA
se_betas <- sumstat$SE
n        <- min(sumstat$N, na.rm = TRUE)

# Extract only genotype dosage columns (assume first column = SNP, rest = genotypes)
#geno_matrix <- as.matrix(dosage[, !colnames(dosage) %in% c("SNP","A1","A2")])
geno_matrix <- as.matrix(dosage[, common_snps])

# Compute LD correlation matrix
R <- cor(geno_matrix, use = "pairwise.complete.obs")
message("✅ Computed LD correlation matrix of dimension: ", nrow(R), "x", ncol(R))

# ---------- Run SuSiE RSS ----------
res_rss <- tryCatch(
  susie_rss(
    bhat = betas,
    shat = se_betas,
    n = n,
    R = R,
    max_iter = susie_iter,
    min_abs_corr = susie_min_abs_cor
  ),
  error = err_handling
)

message("🎉 SuSiE RSS completed successfully.")

# ---------- Extract & Save Results ----------

# Credible sets
cs <- susie_get_cs(res_rss, X = NULL)


if (length(cs$cs) == 0) {
  message("⚠️ No credible sets found for this region.")
  
  } else {
  
    # extract seqid and locus tag from filename (helps concatenation later)
    locuseq <- sub("_sumstat\\.csv$", "", basename(path_sumstat))
    tag_seqid <- sub("_.*$", "", locuseq)
    tag_locus <- sub("^seq\\.\\d+\\.\\d+_", "chr", locuseq)
    
    
    # Build long table of all CS members
    res_list <- lapply(seq_along(cs$cs), function(k) {
      
      idx <- cs$cs[[k]]
      
      # guard against NULL/empty/invalid indices
      if (is.null(idx) || length(idx) == 0) return(NULL)
      idx <- as.integer(idx)
      idx <- idx[is.finite(idx) & idx >= 1 & idx <= nrow(sumstat)]
      if (length(idx) == 0) return(NULL)
      
      data.table(
        seqid = tag_seqid,
        locus = tag_locus,
        cs_id = paste0("CS", k),
        SNP   = sumstat$SNPID[idx],
        CHR   = sumstat$CHR[idx],
        POS   = sumstat$POS[idx],
        BETA  = sumstat$BETA[idx],
        SE    = sumstat$SE[idx],
        PIP   = round(res_rss$pip[idx], 6)
        # A1     = merged$A1_sum[idx],
        # A2     = merged$A2_sum[idx]
      )
    })
    
    cs_summary <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
    
    # save GWAS summary for cs variants
    write.table(
      cs_summary,
      file = out_cs_summary,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
    )
    
    # save full model fitness
    saveRDS(res_rss, file = out_cs_rds)
    
    message("✅ Saved credible set results to: ", out_cs_summary)
}

# full model summary
full_res <- summary(res_rss)

# prepare for merge
cs_details <- full_res$cs %>%
  data.frame() %>%
  mutate(
    seqid = tag_seqid,
    locus = tag_locus,
    cs_id = paste0("CS", cs)
    ) %>%
  select(- cs, - variable)

# list of CS variants
cs_list <- cs_summary %>%
  summarize(cs_snps = paste0(unique(SNP), collapse = ","), .by = "cs_id") %>%
  full_join(cs_details, ., by = "cs_id") %>%
  relocate(seqid, locus, cs_id, cs_log10bf, cs_avg_r2, cs_min_r2, cs_snps)


write.table(cs_list, file = out_cs_list, sep = "\t", row.names = F, quote = F)
message("✅ Saved credible set list to: ", out_cs_list)
message("✅ Analysis done!")
