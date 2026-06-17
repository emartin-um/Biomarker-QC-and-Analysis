#!/usr/bin/env Rscript
# scaffold_qc_run.R — create a fresh, isolated QC run folder under ../../QC_Runs/
#
# WHY: the repo holds the *canonical* pipeline code; each real batch runs in its own folder under
# Alamar_QC_Protocol/QC_Runs/ so runs never overwrite each other or the repo. This script copies the
# CURRENT repo code into a new run folder (freezing the code version that will produce that run's
# outputs — good provenance) and creates empty input_files/ for you to drop the batch data into.
# It copies CODE ONLY (no data, no outputs) so it's fast and can't leak a previous run's data.
#
# USAGE (run from the repo root, i.e. this script's directory):
#   Rscript scaffold_qc_run.R "QC_8plate_10_20_2026"
#   Rscript scaffold_qc_run.R "QC_8plate_10_20_2026" "Primary_QC,Metadata_Merge,Secondary_QC"
#
# Then: drop the batch's NPQ/Sample_QC/Raw_counts into <run>/Primary_QC/input_files/, Knit Primary_QC,
# copy its post-QC NPQ + the metadata into <run>/Metadata_Merge/input_files/, Knit Metadata_Merge.

suppressWarnings(suppressMessages({ library(tools) }))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !nzchar(args[1]))
  stop("Provide a run name, e.g.: Rscript scaffold_qc_run.R \"QC_8plate_10_20_2026\"")

run_name <- args[1]
modules  <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else c("Primary_QC", "Metadata_Merge")
modules  <- trimws(modules)

# ── Resolve paths ─────────────────────────────────────────────────────────────
# Repo root = directory containing this script; QC_Runs sits two levels up (Alamar_QC_Protocol/QC_Runs).
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
REPO <- if (length(this_file)) normalizePath(dirname(this_file)) else normalizePath(getwd())
QC_RUNS <- normalizePath(file.path(REPO, "..", "..", "QC_Runs"), mustWork = FALSE)
RUN_DIR <- file.path(QC_RUNS, run_name)

if (!dir.exists(QC_RUNS)) dir.create(QC_RUNS, recursive = TRUE)
if (dir.exists(RUN_DIR)) stop("Run folder already exists (refusing to overwrite): ", RUN_DIR)

# run tag drives the Metadata_Merge output_dir + the auto-named HTML (e.g. ...Pipeline_8plate_10_20_2026.html)
run_tag <- sub("^QC[_-]?", "", run_name)

# ── Copy CODE ONLY for each module ─────────────────────────────────────────────
# keep code/docs; never copy data, outputs, caches, rendered HTML, RStudio state, or archives.
.keep_pat <- "\\.(Rmd|qmd|R|Rproj|md|yml|yaml)$"
.drop_pat <- paste0("(^|/)(input_files|output_files[^/]*|\\.Rproj\\.user|_cache|.*_files|old|_archive[^/]*|",
                    "_backup[^/]*|example)(/|$)")   # 'example' = synthetic repo demo (e.g. Metadata_Merge/example/), not a real run
copied <- 0L
for (m in modules) {
  src <- file.path(REPO, m)
  if (!dir.exists(src)) { message("! skipping '", m, "' (not found in repo)"); next }
  files <- list.files(src, recursive = TRUE, full.names = FALSE)
  take  <- files[grepl(.keep_pat, files, ignore.case = TRUE) &
                 !grepl(.drop_pat, files, ignore.case = TRUE) &
                 !grepl("\\.(html|RData|rds|csv|xlsx|log|bak)$", files, ignore.case = TRUE)]
  for (f in take) {
    dst <- file.path(RUN_DIR, m, f)
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(src, f), dst, overwrite = TRUE)
    copied <- copied + 1L
  }
  # fresh empty input_files/ at the module root (the pipelines create their own output dirs)
  dir.create(file.path(RUN_DIR, m, "input_files"), recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(RUN_DIR, m, "input_files", ".gitkeep"))
  message("✓ ", m, ": copied ", length(take), " code file(s)")
}

# ── Point Metadata_Merge at this run (output_dir / review_dir) ──────────────────
mm_rmd <- file.path(RUN_DIR, "Metadata_Merge", "Metadata_Merge_Pipeline.Rmd")
if (file.exists(mm_rmd)) {
  L <- readLines(mm_rmd, warn = FALSE)
  L <- sub('^(\\s*output_dir:\\s*).*$', paste0('\\1"output_files_', run_tag, '"   # run-specific (', run_name, ')'), L)
  L <- sub('^(\\s*review_dir:\\s*).*$',  paste0('\\1"review"   # run-specific review dir for ', run_name), L)
  writeLines(L, mm_rmd)
  message("✓ Metadata_Merge: output_dir -> output_files_", run_tag)
}

# ── Provenance stamp ────────────────────────────────────────────────────────────
git_rev <- tryCatch(sub("\\s+$", "", system(paste("git -C", shQuote(REPO), "rev-parse --short HEAD"),
                                            intern = TRUE, ignore.stderr = TRUE)),
                    error = function(e) "unknown")
writeLines(c(
  paste0("Run folder:   ", run_name),
  paste0("Created:      ", as.character(Sys.Date())),
  paste0("Scaffolded from repo: ", REPO),
  paste0("Repo git rev: ", if (length(git_rev)) git_rev else "unknown"),
  paste0("Modules:      ", paste(modules, collapse = ", ")),
  "",
  "Code was copied (frozen) from the repo at the rev above. Re-scaffold to pick up later code changes.",
  "Drop batch data into each module's input_files/ then Knit (Primary_QC first, then Metadata_Merge)."
), file.path(RUN_DIR, "RUN_INFO.txt"))

message("\n✓ Scaffolded ", copied, " code file(s) -> ", RUN_DIR)
message("Next: add data to <run>/Primary_QC/input_files/, Knit Primary_QC, then feed Metadata_Merge.")
