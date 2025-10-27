
# =============================
# SuSiE Fine-mapping Input Loader
# =============================

# Get log path from Snakemake, fallback if missing
log_file <- tryCatch(snakemake@log[[1]], error = function(e) "logs/susieR/default.log")

# Ensure the directory exists
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

# Redirect stdout and stderr to the log file
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")   # redirect stdout
sink(log_con, type = "message")  # redirect messages / stderr

#----------------------------------------#
# ------       Load libraries      ------
#----------------------------------------#

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(data.table)  # For fast I/O
  library(susieR)
  library(pgenlibr) # to load PGEN file
  library(Rfast) #to calculate correlation matrix faster
})

#----------------------------------------#
# -----         User Inputs        ------
#----------------------------------------#

path_sumstat <- snakemake@input[["sumstat"]]
path_pgen <- snakemake@input[["dosage"]]
path_pvar <- gsub(".pgen", ".pvar", path_pgen)
path_psam <- gsub(".pgen", ".psam", path_pgen)

# Load parameters for susieR model
label_chr <- snakemake@params[["chrcol"]]
susie_min_abs_cor <- snakemake@params[["min_abs_corr"]]
susie_iter <- snakemake@params[["iter"]]
#susie_iter <- 1000
#susie_min_abs_cor <- 0.0


# Set TRUE to compute correlation from X, FALSE to load pre-computed LD
compute_ld_from_X <- snakemake@params[["ld_cor"]] # <-- user sets this
path_ld_matrix <- snakemake@input[["ld"]]
path_ld_header <- gsub(".matrix", ".headers", path_ld_matrix)

# outputs
out_data_report <- snakemake@output[["data_report"]]
out_cs_summary <- snakemake@output[["cs_summary"]]
out_cs_list <- snakemake@output[["cs_list"]]
out_cs_rds <- snakemake@output[["cs_rds"]]

#----------------------------------------#
# -------     Helper Functions     ------
#----------------------------------------#

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
check_file(path_pgen)
if (!compute_ld_from_X) { check_file(path_ld_matrix) }


#----------------------------------------#
# ------         Load Data        -------
#----------------------------------------#

# Use fread with explicit arguments to avoid surprises
sumstat <- tryCatch({
  fread(path_sumstat, header = TRUE, sep = "\t", data.table = FALSE)
}, error = function(e) {
  stop("❌ Failed to read sumstat file: ", e$message)
})

# number of SNPs in GWAS results subset
n_snp_sumstat <- nrow(sumstat)

# Read psam and pvar
psam_df <- read.delim(path_psam, header = TRUE, comment.char = "")
pvar_df <- read.delim(path_pvar, header = TRUE, comment.char = "")

# Read pgen
pgen <- tryCatch({
  # Read pgen
  #pvar <- pgenlibr::NewPvar(path_pvar)
  NewPgen(path_pgen) #, pvar=pvar
  }, error = function(e) {
    stop("❌ Failed to read dosage file: ", e$message)
})

#----------------------------------------#
# --------       Basic QC         -------
#----------------------------------------#
# rename column name
colnames(sumstat)[which(names(sumstat) == label_chr)] <- "CHR"

# Check mandatory columns in summary stats
required_sumstat_cols <- c("SNPID", "CHR", "POS", "EA", "NEA", "BETA", "SE", "MLOG10P")
missing_cols <- setdiff(required_sumstat_cols, colnames(sumstat))

if (length(missing_cols) > 0) {
  stop("❌ Missing required columns in sumstat file: ", paste(missing_cols, collapse = ", "))
}

#-------------# 
# Check the number of variants and samples
n_variants <- pgenlibr::GetVariantCt(pgen)
n_samples  <- pgenlibr::GetRawSampleCt(pgen)

# Extract dosages for all of variants
dosage <- pgenlibr::ReadList(pgen, 1:n_variants, meanimpute = FALSE)

# Add variant IDs as column names
colnames(dosage) <- pvar_df$ID

# Add sample IDs as row names
rownames(dosage) <- psam_df$IID


message("✅ Summary stats and dosage files loaded successfully.")

#----------------------------------------#
# -------     Variant Matching     ------
#----------------------------------------#

# to avoid allele mismatch
common_snps <- intersect(sumstat$SNPID, pvar_df$ID)
n_common_snps <- length(common_snps)

if (n_common_snps == 0) {
  stop("❌ No overlapping SNPs between sumstat and dosage files.")
}

message("✅ Overlapping SNPs found: ", n_common_snps)

# Optional: subset both datasets to common SNPs
sumstat <- sumstat[sumstat$SNPID %in% common_snps, ]
X <- dosage[, common_snps] %>% as.matrix()

message("✅ Subsetted to common SNPs. Ready for SuSiE.")

#----------------------------------------#
# ---------- Allele Alignment ----------
#----------------------------------------#

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

#----------------------------------------#
# ------      Reporting Counts      -----
#----------------------------------------#

#extract seqid_locus compound
locuseq <- sub("_sumstat\\.csv$", "", basename(path_sumstat))

# reporting numbers of input data 
data_counts <- data.frame(
  "locuseq"       = locuseq,
  "n_snp_pgen"    = n_variants,
  "n_snp_gwas"    = n_snp_sumstat,
  "n_snp_shared"  = n_common_snps,
  "n_sample_pgen" = n_samples
)

# saving the report
write.table(data_counts, file = out_data_report, sep = "\t", row.names = F, quote = F)

message("✅ Saved  input characteristics  to : ", out_data_report)


#----------------------------------------#
# -----  Load or Compute LD matrix  -----
#----------------------------------------#

if (compute_ld_from_X) {
  message("📈 Computing LD correlation matrix from genotype matrix X ...")
  R <- cor(X, use = "pairwise")
  } else {
    message("📥 Loading precomputed LD matrix from PLINK2 output ...")
    ld_headers <- fread(path_ld_header, header = FALSE, col.names = "SNP")
    R <- fread(path_ld_matrix) %>% as.matrix()
    rownames(R) <- colnames(R) <- ld_headers$SNP
    }

message("✅ LD matrix of dimension: ", nrow(R), "x", ncol(R))

# CHECK SYMMETRY: 
if (!isSymmetric(R, tol = 1e-8)) {
  stop("❌ The LD matrix is not symmetric. Please check your input.")
}

# Check Positive semi-definiteness: Cholesky factorization 
# of a real symmetric positive-definite square matrix
positive_semi_definite <- TRUE
tryCatch({
  invisible(chol(R))   # will fail if not positive-definit
  }, error = function(e) {
    positive_semi_definite <- FALSE
})

if (!positive_semi_definite) {
  stop("❌ The LD matrix is not positive semi-definite. SuSiE requires PSD LD matrix.")
} else {
  message("✅ The LD matrix is symmetric and positive semi-definite.")
}


#----------------------------------------#
# -----  Prepare Inputs for SuSiE   -----
#----------------------------------------#

betas    <- sumstat$BETA
se_betas <- sumstat$SE
n        <- min(sumstat$N, na.rm = TRUE)

#----------------------------------------#
# -------      Run SuSiE RSS      -------
#----------------------------------------#

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

#----------------------------------------#
# -----    Extract & Save Results  ------
#----------------------------------------#

# full model summary
full_res <- summary(res_rss)

# Credible sets
#cs <- susie_get_cs(res_rss, X = NULL) # issue #257: generates more credible sets as it applies NO impurity filter
cs <- full_res$cs

if (length(cs$cs) == 0) {
  message("⚠️ No credible sets found for this region.")
  
  } else {
  
    # extract seqid and locus tag from filename (helps concatenation later)
    tag_seqid <- sub("_.*$", "", locuseq)
    tag_locus <- sub("^seq\\.\\d+\\.\\d+_", "chr", locuseq)
    
    
    # Build long table of all CS members
    res_list <- lapply(seq_along(cs$cs), function(k) {
      
      # CS index
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
    message("✅ Saved SuSiE full summary to: ", out_cs_rds)
}


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
  summarize(
    cs_snps = paste0(unique(SNP), collapse = ","), 
    .by = c("seqid", "locus", "cs_id") # keep seqid and locus in CS list
    ) %>%
  full_join(cs_details, ., join_by(seqid, locus, cs_id)) %>%
  relocate(seqid, locus, cs_id, cs_log10bf, cs_avg_r2, cs_min_r2, cs_snps)


write.table(cs_list, file = out_cs_list, sep = "\t", row.names = F, quote = F)

message("✅ Saved credible set list to: ", out_cs_list)
message("✅ Analysis done!")


#-------------# 
# Reset sinks at the end
sink(type = "message")
sink(type = "output")
close(log_con)
