
# =============================
# SuSiE Fine-mapping Log
# =============================

# Get log path from Snakemake, fallback if missing
log_file <- tryCatch(snakemake@log[[1]], error = function(e) "logs/susieR/default.log")

# Ensure the directory exists
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

# Redirect stdout and stderr to the log file
log_con <- file(log_file, open = "wt")

# Safety reset in case of prior crash
tryCatch({
  while (sink.number()  > 0) sink()
}, error = function(e) invisible(NULL))

sink(log_con, type = "output")   # redirect stdout
sink(log_con, type = "message")  # redirect messages / stderr

# Reset sinks at the end
on.exit({
  sink(type = "message")   # restore message sink first
  sink(type = "output")    # then output sink
  close(log_con)
}, add = TRUE)


# Suppress warning for kriging plot
options(lifecycle_verbosity = "quiet")  

#-------------#

start_time <- Sys.time()
start_time

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
  library(coloc)
  library(ggplot2)
  library(ggplotify)
  library(ggcorrplot)
  library(ggpubr)
})


set.seed(777)

#----------------------------------------#
# -----         User Inputs        ------
#----------------------------------------#

path_sumstat <- snakemake@input[["sumstat"]]
path_pgen <- snakemake@input[["pgen"]]
path_pvar <- snakemake@input[["pvar"]]
path_psam <- snakemake@input[["psam"]]

# Load parameters for susieR model
label_chr <- snakemake@params[["chrcol"]]
susie_min_abs_cor <- snakemake@params[["min_abs_corr"]]
susie_iter <- snakemake@params[["iter"]]
susie_L <- snakemake@params[["L"]]
susie_est_resvar <- snakemake@params[["est_res_var"]]
study_id <- snakemake@params[["study"]]
n_gwas <- snakemake@params[["n_gwas"]]


# Set TRUE to compute correlation from X, FALSE to load pre-computed LD
compute_ld_from_X <- snakemake@params[["ld_cor"]] # <-- user sets this
path_ld_matrix <- snakemake@input[["ld"]]
path_ld_header <- snakemake@input[["ld_snps"]]

# outputs
out_data_report <- snakemake@output[["data_report"]]
out_cs_summary <- snakemake@output[["cs_summary"]]
out_cs_list <- snakemake@output[["cs_list"]]
out_cs_rds <- snakemake@output[["cs_rds"]]
out_cs_annot <- snakemake@output[["cs_annot"]]
out_kriging <- snakemake@output[["kriging"]]


# Extract tags from filename helps concatenation later
locuseq <- sub("_sumstat\\.csv$", "", basename(path_sumstat))
tag_seqid <- sub("_.*$", "", locuseq)
tag_locus <- sub("^seq\\.\\d+\\.\\d+_", "chr", locuseq)


#----------------------------------------#
# -------     Helper Functions     ------
#----------------------------------------#

susie_warnings <- character()

wrn_handling <- function(w) {
  cleaned_warn <- gsub("\033\\[[0-9;]*m|\\n", "", conditionMessage(w))
  susie_warnings <<- c(susie_warnings, cleaned_warn)
  message("⚠️ ", conditionMessage(w))
  invokeRestart("muffleWarning")
}

err_handling <- function(e) { stop("❌ SuSiE failed: ", e$message) }


check_file <- function(path, min_size = 1e2) {
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

if (compute_ld_from_X){
  check_file(path_pgen)
  } else {
    check_file(path_ld_matrix)
  }


#----------------------------------------#
# ------         Load Data        -------
#----------------------------------------#


# Use fread with explicit arguments to avoid surprises
if (study_id == "interval") {
  
  headers <- c("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "N", "BETA", "SE", "MLOG10P", "CHISQ")
  
  sumstat <- tryCatch({
    fread(path_sumstat, header = FALSE, col.names = headers, sep = "\t", data.table = FALSE)
    }, error = function(e) {
      stop("❌ Failed to read sumstat file: ", e$message)
      })
  
  } else if (study_id == "believe") {
    
    headers <- c("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "BETA", "SE", "P", "MLOG10P", "Z")
    
    sumstat <- tryCatch({
      fread(path_sumstat, header = FALSE, col.names = headers, sep = "\t", data.table = FALSE)
      }, error = function(e) {
        stop("❌ Failed to read sumstat file: ", e$message)
        })
    
    # define sample size manually
    sumstat$N <- n_gwas
    
    } else if (study_id == "meta") {
      
      sumstat <- tryCatch({
        fread(path_sumstat, header = TRUE, sep = "\t", data.table = FALSE)
        }, error = function(e) {
          stop("❌ Failed to read sumstat file: ", e$message)
          })
      
      } else {
        stop("❌ The study name is incorrect. Please use any of interval/believe/meta in config file.")
}

message("✅ Summary stats and variant files loaded successfully.")


#----------------------------------------#
# --------       Basic QC         -------
#----------------------------------------#

# rename column name
if (study_id == "meta") {
  colnames(sumstat)[which(names(sumstat) == label_chr)] <- "CHR"
}

# Check mandatory columns in summary stats
required_sumstat_cols <- c("SNPID", "CHR", "POS", "EA", "NEA", "BETA", "SE", "MLOG10P", "N")
missing_cols <- setdiff(required_sumstat_cols, colnames(sumstat))

if (length(missing_cols) > 0) {
  stop("❌ Missing required columns in sumstat file: ", paste(missing_cols, collapse = ", "))
}


#----------------------------------------#
# -------     Variant Matching     ------
#----------------------------------------#

if (compute_ld_from_X) {
  
  ld_size <- NA_character_
  
  # Read pgen
  pgen <- tryCatch({
    #pvar <- pgenlibr::NewPvar(path_pvar)
    pgenlibr::NewPgen(path_pgen) #, pvar=pvar
    }, error = function(e) {
      stop("❌ Failed to read dosage file: ", e$message) 
    })
  
  # Read psam & pvar
  psam_df <- read.delim(path_psam, header = TRUE, comment.char = "")
  pvar_df <- read.delim(path_pvar, header = TRUE, comment.char = "")
  
  # Cross variants in GWAS & PGEN/LD to avoid allele mismatch
  common_snps  <- intersect(sumstat$SNPID, pvar_df$ID)
  
  } else {
    
    # Report LD file size along with data counts
    ld_size <- round(file.size(path_ld_matrix)/(1024*1024), 2)
  
    # Read SNPs in pre-computed LD file
    ld_header <- fread(path_ld_header, header = FALSE, col.names = "SNP")
    
    R <- tryCatch({
      fread(path_ld_matrix, header = FALSE, sep = "\t", data.table = FALSE)
    }, error = function(e) {
      stop("❌ Failed to read LD file: ", e$message)
    })
    
    R <- as.matrix(R)
    rownames(R) <- colnames(R) <- ld_header$SNP
    
    # Cross variants to avoid allele mismatch
    common_snps  <- intersect(sumstat$SNPID, ld_header$SNP)

}

#-------------#

n_snp_common <- length(common_snps)

if (n_snp_common == 0) {
  stop("❌ No overlapping SNPs between sumstat and dosage files.")
}


#----------------------------------------#
# -----  Load or Compute LD matrix  -----
#----------------------------------------#

if (compute_ld_from_X) {
  
  message("📈 Computing LD correlation matrix from genotype matrix X ...")
  
  # Count number of variants & samples in raw data
  n_variants <- nrow(pvar_df) # Or → pgenlibr::GetVariantCt(pgen)
  n_samples  <- nrow(psam_df) # Or → pgenlibr::GetRawSampleCt(pgen)
  
  # Extract dosages for all of variants
  dosage <- pgenlibr::ReadList(pgen, 1:n_variants, meanimpute = FALSE)
  
  # Add variant IDs as column names
  colnames(dosage) <- pvar_df$ID
  
  # Add sample IDs as row names
  rownames(dosage) <- psam_df$IID
  
  # Define genotype matrix using shared variants
  #   - Columns: variants dosage levels
  #   - Rows: all individuals in PGEN file
  X <- dosage[, common_snps] %>% as.matrix()
  
  # Compute LD (r) correlation matrix
  R <- cor(X, use = "pairwise")
  
  } else{
    
    message("📥 Loading precomputed LD matrix from PLINK2 output ...")
    
    n_variants <- NA_character_
    n_samples  <- NA_character_
    
    # Select only variants in region sumstat
    R <- R[common_snps, common_snps]

}

#-------------#

# Number of variants in GWAS
n_snp_sumstat <- nrow(sumstat)

# Subset both datasets to common SNPs
sumstat <- sumstat[sumstat$SNPID %in% common_snps, ]

message("✅ Overlapping SNPs found: ", n_snp_common)
message("✅ GWAS & LD were limited to common SNPs.")
message("✅ LD matrix of dimension: ", nrow(R), "x", ncol(R))


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
# -----        Count Variants       -----
#----------------------------------------#

# Define variant type per row
sumstat <- sumstat %>%
  mutate(
    is_indel = nchar(EA) != 1 | nchar(NEA) != 1,
    is_snp   = nchar(EA) == 1 & nchar(NEA) == 1
  )

# Count indels and SNPs
n_indels <- sum(sumstat$is_indel)
n_snps   <- sum(sumstat$is_snp)

# Site-level allele counts
site_counts <- sumstat %>%
  group_by(CHR, POS) %>%
  summarise(
    n_alleles = n_distinct(c(EA, NEA)),
    .groups = "drop"
  )

# Count bi-allelic and multi-allelic sites
n_biallelic    <- sum(site_counts$n_alleles == 2)
n_multiallelic <- sum(site_counts$n_alleles > 2)


#----------------------------------------#
# -----      Inspect LD matrix      -----
#----------------------------------------#

# Inspect LD matrix symmetry 
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
z_scores <- betas / se_betas
n        <- min(sumstat$N, na.rm = TRUE)


#----------------------------------------#
# ----   Quantify LD Misalignment    ----
#----------------------------------------#

warntxt <- NA_character_

# Capture warning while estimating lambda
withCallingHandlers(
  {
    lambda <- estimate_s_rss(z = z_scores, R = R) # n = n adjusts for sample size
  },
  message = function(m) {                 # Let it continue to the sink
    warntxt <<- gsub("\033\\[[0-9;]*m|\\n", "", conditionMessage(m)) # omit HTML coloring warning
    message("⚠️ ", conditionMessage(m))   # re-emit original (with color) → goes to log
    invokeRestart("muffleMessage")        # suppress duplicate print
  }
)

message("✅ Estimated λ measuring LD mismatch: ", signif(lambda, 4))

#-------------#

# Conditional z-score diagnostics
withCallingHandlers(
  {
    condz <- kriging_rss(z = z_scores, R = R) # n = n adjusts for sample size
  },
  message = function(m) invokeRestart("muffleMessage")
  )

# Add lambda annotation to diagnostic plot
plt_kriging <- condz$plot +
  annotate(
    "text",
    label = paste0("λ = ", signif(lambda, 4)),
    hjust = "left",
    vjust = "top",
    size = 4,
    x = -Inf,
    y = Inf
  ) +
  labs(
    title = "LD matrix consistency with Z-scores"
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 10, face = 2, hjust = 0.5),
    axis.title = element_text(size = 12),
    panel.background = element_blank()
  )


#----------------------------------------#
# ------      Reporting Counts      -----
#----------------------------------------#

# Report input data counts
data_counts <- data.frame(
  "seqid"          = tag_seqid,
  "locus"          = tag_locus,
  "nsample_pgen"   = n_samples,
  #"nsample_gwas"  = n_gwas,
  "nvar_pgen"      = n_variants,
  "nvar_gwas"      = n_snp_sumstat,
  "nvar_shared"    = n_snp_common,
  "nsnps_shared"   = n_snps,
  "nindels_shared" = n_indels,
  "nbi_allelic"    = n_biallelic,
  "nmulti_allelic" = n_multiallelic,
  "run_time_min"   = NA_character_,
  "ld_from_X"      = compute_ld_from_X,
  "ld_size_mg"     = ld_size,
  "lambda"         = lambda,
  "lambda_warning" = warntxt,
  "susie_warning"  = NA_character_
)


#----------------------------------------#
# -------      Run SuSiE RSS      -------
#----------------------------------------#

message("▶ Running SuSiE...")

res_rss <- tryCatch(
  withCallingHandlers(
    {
      susie_rss(
        z = z_scores,
        z_method = "score", # "score" assumes z (from mixed model) are on PVE scale, so no adjustment is applied.
        n = n,
        R = R,
        L = susie_L,
        max_iter = susie_iter,
        min_abs_corr = susie_min_abs_cor,
        estimate_residual_variance = susie_est_resvar # TRUE if using in-sample LD
      )
    }, warning = wrn_handling
  ), error = err_handling
)

message("🎉 SuSiE completed.")

#----------------------------------------#
# -----    Extract & Save Results  ------
#----------------------------------------#

# full model summary
full_res <- summary(res_rss)

# Credible sets
#cs <- susie_get_cs(res_rss, X = NULL) # issue #257: generates more credible sets as it applies NO impurity filter

cs   <- full_res$cs    # containing CS impurity indices
vars <- full_res$vars  # containing CS Posterior Inclusion Probabilities

# Handle regions with no credible sets
if (is.null(cs$cs) || length(cs$cs) == 0) {
  message("⚠️ No credible sets found for region ", locuseq)
  
  # create a fake GWAS summary
  cs_summary <- data.table(
    seqid = tag_seqid,
    locus = tag_locus,
    cs_id = "no_credible",
    PIP   = NA_character_
  )
  
  # create a fake list of CS variants
  cs_list <- data.table(
    seqid = tag_seqid,
    locus = tag_locus,
    cs_id = "no_credible",
    cs_log10bf = NA_character_,
    cs_avg_r2 = NA_character_,
    cs_min_r2 = NA_character_,
    ncs = NA_character_,
    cs_snps = NA_character_
  )
  
  plt_cs <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "No credible set identified",
      size = 5,
      hjust = 0.5,
      vjust = 0.5
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      title = "Posterior Inclusion Probablity of variants"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  } else {
    
    # list of the entire SNPs with PIP
    snps_pip <- vars %>%
      transmute(
        cs_id = cs,
        SNPID = sumstat$SNPID[variable],
        PIP = variable_prob
      )
    
    # subset of GWAS results for CS variants
    cs_summary <- sumstat %>%
      left_join(snps_pip, by = "SNPID") %>%
      left_join(cs[1:4], join_by(cs_id == cs)) %>% # append CS characteristics to sumstat
      filter(cs_id > 0)
    
    cs_sum_plot <- cs_summary %>% mutate(
      cs_id = factor(
        cs_id,
        levels = sort(unique(cs_id)),
        labels = paste0("L", sort(unique(cs_id))),
        ordered = TRUE
      ))
    
    # list of CS variants
    cs_list <- cs_summary %>%
      summarize(
        cs_snps = paste(SNPID, collapse = ","),
        .by = cs_id
        ) %>%
      full_join(
        cs,
        join_by(cs_id == cs)) %>% # add impurity indices
      mutate(
        seqid = tag_seqid, # store seqid and locus in CS list
        locus = tag_locus,
        ncs = str_count(variable, ",") + 1  # number of variants in each set
      ) %>% # only removed 'variable', indicating CS indices
      select(seqid, locus, cs_id, cs_log10bf, cs_avg_r2, cs_min_r2, ncs, cs_snps)
    
    
    # plot credible sets
    plt_cs <- ggplotify::as.ggplot(function() {
      
      par(
        mar = c(7,5,3,3),
        bg = "white",
        fg = "black"
        )
      
      susie_plot(
        res_rss,
        y = "PIP",
        b = betas,
        xlab = "Variants",
        add_bar = FALSE,
        add_legend = TRUE,
        cex.axis = .75
        #main = paste("SeqID:", tag_seqid, "\nRegion:", tag_locus)
      )
    })
    
    plt_cs <- plt_cs +
      labs(
        title = "Posterior Inclusion Probablity of variants"
      ) +
      theme(
        plot.title = element_text(size = 10, face = 2, hjust = 0.5),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA)
      )
    
}


#----------------------------------------#
# ------       Visualizations      ------
#----------------------------------------#

# Regional plot
plt_lz <- sumstat %>%
  ggplot(aes(x = POS, y = MLOG10P)) +
  geom_point(size = 3, color = "#FDC700", shape = 16, alpha = .7) +
  geom_point(data = cs_sum_plot, aes(color = cs_id), size = 4, shape = 21, stroke = 2) +
  scale_x_continuous(labels = function(x) round(x/1e6, 2)) +
  labs(
    title = paste("Region:", locuseq),
    x = "Genomic Position (Mb)"
  ) +
  theme_light() +
  theme(
    legend.position = c(.94, .85),
    legend.background = element_blank(),
    plot.title = element_text(size = 10, face = 2, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 9),
    axis.ticks.length = unit(2.5, "mm")
  )


# CS LD Correlation
# Flag one credible sets
if (length(cs$cs) >= 2) {
  
  # Calculate CS correlation
  cs_ld <- susieR::get_cs_correlation(res_rss, Xcorr = R)
  
  # Numerical ordering: L1, L2, L3, ...
  cs_order <- order(
    as.numeric(sub("^L", "", rownames(cs_ld)))
  )
  
  cs_ld <- cs_ld[cs_order, cs_order]
  
  # Correlation plot
  plt_cs_ld <- ggcorrplot::ggcorrplot(
    cs_ld,
    method = "square",
    type = "lower",
    lab = TRUE,
    lab_size = 3,
    tl.cex = 10, # axis text
    digits = 2,
    colors = c("#3B4CC0", "white", "#B40426"),
    outline.color = "white",
    ggtheme = ggplot2::theme_minimal()
  ) +
    labs(
      title = "LD correlation between credible sets",
      x = NULL,
      y = NULL
    ) +
    theme(
      plot.title = element_text(size = 10, face = 2, hjust = 0.5),
      panel.background = element_blank()
    )
  
  } else {
  
  # Only one credible set → no pairwise LD correlation exists
  plt_cs_ld <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "Only one credible set\nNo pairwise LD correlation",
      size = 5,
      hjust = 0.5,
      vjust = 0.5
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      title = "LD correlation between credible sets"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14)
    )
}


# Combine plots
plt_joint <- ggpubr::ggarrange(
  plt_lz,
  plt_kriging,
  plt_cs,
  plt_cs_ld,
  nrow = 2,
  ncol = 2,
  #align = "hv",
  heights = c(1, 1)
)


# Store multi-panel plot
ggsave(filename = out_kriging, plt_joint, width = 9, height = 7.75, dpi = 200)
message("✅ Saved susie visualizations to: ", out_kriging)


#-------------#
# save GWAS summary for cs variants
write.table(cs_summary, file = out_cs_summary, sep = "\t", row.names = F, quote = FALSE)
message("✅ Saved credible set results to: ", out_cs_summary)

# save full model fitness
saveRDS(res_rss, file = out_cs_rds)
message("✅ Saved SuSiE full summary to: ", out_cs_rds)

# Annotate and save full model fitness with LD matrix for coloc
res_rss_annot <- coloc::annotate_susie(res_rss, sumstat$SNPID, R)

saveRDS(res_rss_annot, file = out_cs_annot)
message("✅ Saved LD-annotated SuSiE full summary for coloc: ", out_cs_annot)

write.table(cs_list, file = out_cs_list, sep = "\t", row.names = F, quote = F)
message("✅ Saved credible set list to: ", out_cs_list)
message("✅ Analysis done!")


#-------------#
# Report run time
end_time <- Sys.time()
end_time

elapsed_time <- round(as.numeric(end_time - start_time, units="mins"), 3)

message("Run time: ", elapsed_time, " minutes\n")

# Append runtime to report
data_counts$run_time_min <- elapsed_time

data_counts$susie_warning <- if (length(susie_warnings) > 0) {
  paste(susie_warnings, collapse = " | ")
  } else {
  NA_character_
}

# saving the report
write.table(
  data_counts,
  file = out_data_report,
  sep = "\t",
  row.names = F,
  quote = T
  )

message("✅ Saved report to: ", out_data_report)

