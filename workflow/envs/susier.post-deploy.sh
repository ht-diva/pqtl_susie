#!/usr/bin/env bash
set -euo pipefail

# Prevent cluster-wide R settings from overriding the Conda R environment.
unset R_HOME
unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

CONDA_R_LIBRARY="${CONDA_PREFIX}/lib/R/library"
mkdir -p "${CONDA_R_LIBRARY}"

echo "Rscript: $(command -v Rscript)"
echo "CONDA_PREFIX: ${CONDA_PREFIX}"
echo "Conda R library: ${CONDA_R_LIBRARY}"

Rscript --vanilla -e '
  conda_lib <- file.path(
    Sys.getenv("CONDA_PREFIX"),
    "lib",
    "R",
    "library"
  )

  if (!dir.exists(conda_lib)) {
    dir.create(conda_lib, recursive = TRUE)
  }

  .libPaths(c(conda_lib, .libPaths()))

  message("R version: ", R.version.string)
  message("R home: ", R.home())
  message("R executable: ", file.path(R.home("bin"), "R"))
  message("Installation library: ", conda_lib)
  message("Library search paths:")
  print(.libPaths())

  if (!file.access(conda_lib, 2) == 0) {
    stop("Conda R library is not writable: ", conda_lib)
  }

  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "false")

  remotes::install_github(
    "RcppCore/RcppParallel",
    lib = conda_lib,
    upgrade = "never",
    force = TRUE
  )

  remotes::install_github(
    "stephenslab/susieR",
    lib = conda_lib,
    upgrade = "never",
    force = TRUE
  )

  message(
    "RcppParallel version: ",
    packageVersion("RcppParallel", lib.loc = conda_lib)
  )
  message(
    "susieR version: ",
    packageVersion("susieR", lib.loc = conda_lib)
  )

  stopifnot(
    packageVersion("RcppParallel", lib.loc = conda_lib) >= "5.1.10",
    packageVersion("susieR", lib.loc = conda_lib) >= "0.16.6"
  )
'