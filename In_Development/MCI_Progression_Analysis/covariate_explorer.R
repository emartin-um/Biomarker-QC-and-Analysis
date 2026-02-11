################################################################################
# covariate_explorer.R
# Standalone utility for exploring, configuring, and applying covariate groupings
#
# Usage: source("covariate_explorer.R") in any R project
# Dependency: tidyverse
################################################################################

library(tidyverse)

################################################################################
# SECTION A: EXPLORATION FUNCTIONS
################################################################################

#' Explore frequency distribution of a single covariate
#'
#' @param data A dataframe
#' @param var_name Character: column name to explore
#' @return Invisible tibble with level, n, pct columns
explore_covariate <- function(data, var_name) {
  stopifnot(var_name %in% names(data))

  counts <- data %>%
    count(!!sym(var_name), name = "n") %>%
    arrange(desc(n)) %>%
    mutate(pct = round(n / sum(n) * 100, 1))

  message(sprintf("\n=== %s Distribution (N=%d) ===", var_name, nrow(data)))
  print(as.data.frame(counts), row.names = FALSE)

  invisible(counts)
}

#' Two-way cross-tabulation of two covariates
#'
#' @param data A dataframe
#' @param row_var Character: column for rows
#' @param col_var Character: column for columns
#' @return Invisible tibble in wide format with row totals
explore_cross_tab <- function(data, row_var, col_var) {
  stopifnot(row_var %in% names(data), col_var %in% names(data))

  ct <- data %>%
    count(!!sym(row_var), !!sym(col_var)) %>%
    pivot_wider(names_from = all_of(col_var), values_from = n, values_fill = 0) %>%
    mutate(Total = rowSums(across(where(is.numeric)))) %>%
    arrange(desc(Total))

  message(sprintf("\n=== %s x %s Cross-Tabulation ===", row_var, col_var))
  print(as.data.frame(ct), row.names = FALSE)

  invisible(ct)
}

#' Explore all specified covariates: individual frequencies + pairwise cross-tabs
#'
#' @param data A dataframe
#' @param covariate_cols Character vector: column names to explore
#' @return Invisible list of all tables produced
explore_all_covariates <- function(data, covariate_cols) {
  results <- list()

  # Individual frequencies
  message("\n##########################################################")
  message("# INDIVIDUAL COVARIATE DISTRIBUTIONS")
  message("##########################################################")
  for (v in covariate_cols) {
    if (v %in% names(data)) {
      results[[v]] <- explore_covariate(data, v)
    } else {
      message(sprintf("\nWARNING: Column '%s' not found in data, skipping.", v))
    }
  }

  # Pairwise cross-tabs
  valid_cols <- covariate_cols[covariate_cols %in% names(data)]
  if (length(valid_cols) >= 2) {
    message("\n##########################################################")
    message("# PAIRWISE CROSS-TABULATIONS")
    message("##########################################################")
    pairs <- combn(valid_cols, 2, simplify = FALSE)
    for (pair in pairs) {
      key <- paste(pair, collapse = "_x_")
      results[[key]] <- explore_cross_tab(data, pair[1], pair[2])
    }
  }

  invisible(results)
}

#' Show all unique combinations of specified columns with counts
#'
#' @param data A dataframe
#' @param ... Column names (unquoted or character)
#' @return Invisible tibble of combinations sorted by descending count
explore_combinations <- function(data, ...) {
  cols <- c(...)
  stopifnot(all(cols %in% names(data)))

  combos <- data %>%
    count(across(all_of(cols)), name = "n") %>%
    arrange(desc(n)) %>%
    mutate(pct = round(n / sum(n) * 100, 1))

  col_label <- paste(cols, collapse = " x ")
  message(sprintf("\n=== All Combinations: %s (N=%d, %d unique) ===",
                  col_label, nrow(data), nrow(combos)))
  print(as.data.frame(combos), row.names = FALSE)

  invisible(combos)
}

################################################################################
# SECTION B: CONFIGURATION HELPERS
################################################################################

#' Create a category config (for recoding a single column)
#'
#' @param source_col Character: name of column to recode
#' @param target_col Character: name of new column to create
#' @param keep Character vector: levels to keep as-is
#' @param collapse Named list: target_label = c(source_levels)
#' @param exclude Character vector: levels to drop entirely
#' @return A category config list
make_category_config <- function(source_col, target_col,
                                  keep = character(0),
                                  collapse = list(),
                                  exclude = character(0)) {
  config <- list(
    type = "category",
    source_col = source_col,
    target_col = target_col,
    keep = keep,
    collapse = collapse,
    exclude = exclude
  )
  class(config) <- c("covariate_config", "category_config", "list")
  config
}

#' Create a grouping config (for combining multiple columns into one)
#'
#' @param target_col Character: name of new column to create
#' @param rules Named list of rules. Each rule is a named list of column=value pairs.
#'   Rules are evaluated in order; first match wins.
#' @param other_label Character: label for rows matching no rule
#' @param exclude_other Logical: if TRUE, drop rows that match no rule
#' @return A grouping config list
make_grouping_config <- function(target_col, rules,
                                  other_label = "Other",
                                  exclude_other = FALSE) {
  config <- list(
    type = "grouping",
    target_col = target_col,
    rules = rules,
    other_label = other_label,
    exclude_other = exclude_other
  )
  class(config) <- c("covariate_config", "grouping_config", "list")
  config
}

#' Pretty-print a covariate config
#'
#' @param config A config list (category or grouping)
print_config <- function(config) {
  if (isTRUE(config$type == "category")) {
    message(sprintf("\n=== Category Config: %s -> %s ===",
                    config$source_col, config$target_col))
    if (length(config$keep) > 0) {
      message(sprintf("  KEEP as-is: %s", paste(config$keep, collapse = ", ")))
    }
    if (length(config$collapse) > 0) {
      for (target in names(config$collapse)) {
        message(sprintf("  COLLAPSE to '%s': %s",
                        target, paste(config$collapse[[target]], collapse = ", ")))
      }
    }
    if (length(config$exclude) > 0) {
      message(sprintf("  EXCLUDE (drop rows): %s", paste(config$exclude, collapse = ", ")))
    }
  } else if (isTRUE(config$type == "grouping")) {
    message(sprintf("\n=== Grouping Config -> %s ===", config$target_col))
    for (rule_name in names(config$rules)) {
      conditions <- config$rules[[rule_name]]
      cond_str <- paste(names(conditions), "=", unlist(conditions), collapse = " AND ")
      message(sprintf("  '%s': where %s", rule_name, cond_str))
    }
    message(sprintf("  Unmatched rows -> '%s' (exclude_other = %s)",
                    config$other_label, config$exclude_other))
  } else {
    message("Unknown config type")
    str(config)
  }
}

################################################################################
# SECTION C: APPLICATION FUNCTIONS
################################################################################

#' Apply a category config to recode a single column
#'
#' @param data A dataframe
#' @param config A category config list
#' @return Modified dataframe with new target column; excluded rows removed
apply_category_config <- function(data, config) {
  stopifnot(config$source_col %in% names(data))

  src <- config$source_col
  tgt <- config$target_col

  # Check for unaccounted levels
  observed <- unique(data[[src]])
  observed <- observed[!is.na(observed)]
  accounted <- c(config$keep, unlist(config$collapse), config$exclude)
  unaccounted <- setdiff(observed, accounted)

  if (length(unaccounted) > 0) {
    warning(sprintf("Unaccounted %s levels (will be dropped): %s",
                    src, paste(unaccounted, collapse = ", ")))
  }

  # Build mapping: source_value -> target_value
  mapping <- character(0)
  for (level in config$keep) {
    mapping[level] <- level
  }
  for (target_label in names(config$collapse)) {
    for (source_level in config$collapse[[target_label]]) {
      mapping[source_level] <- target_label
    }
  }

  # Apply
  n_before <- nrow(data)
  data[[tgt]] <- mapping[data[[src]]]

  # Remove excluded and unaccounted (NA in target)
  data <- data[!is.na(data[[tgt]]), ]
  n_after <- nrow(data)

  message(sprintf("\n--- apply_category_config: %s -> %s ---", src, tgt))
  message(sprintf("  Rows: %d -> %d (removed %d)", n_before, n_after, n_before - n_after))
  message(sprintf("  %s distribution:", tgt))
  ct <- sort(table(data[[tgt]], useNA = "ifany"), decreasing = TRUE)
  for (i in seq_along(ct)) {
    message(sprintf("    %-30s %d", names(ct)[i], ct[i]))
  }

  data
}

#' Apply a grouping config to create a new column from multi-column rules
#'
#' @param data A dataframe
#' @param config A grouping config list
#' @return Modified dataframe with new target column; optionally filtered
apply_grouping_config <- function(data, config) {
  tgt <- config$target_col

  # Validate that all rule columns exist
  all_rule_cols <- unique(unlist(lapply(config$rules, names)))
  missing_cols <- setdiff(all_rule_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Columns referenced in rules but not in data: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  # Initialize to other_label
  data[[tgt]] <- config$other_label

  # Apply rules in REVERSE order so first rule in list has highest priority
  for (i in rev(seq_along(config$rules))) {
    rule_name <- names(config$rules)[i]
    rule <- config$rules[[i]]

    mask <- rep(TRUE, nrow(data))
    for (col_name in names(rule)) {
      col_match <- data[[col_name]] == rule[[col_name]]
      col_match[is.na(col_match)] <- FALSE
      mask <- mask & col_match
    }

    data[[tgt]][mask] <- rule_name
  }

  # Optionally exclude unmatched
  n_before <- nrow(data)
  if (config$exclude_other) {
    data <- data[data[[tgt]] != config$other_label, ]
  }
  n_after <- nrow(data)

  message(sprintf("\n--- apply_grouping_config -> %s ---", tgt))
  if (n_before != n_after) {
    message(sprintf("  Rows: %d -> %d (removed %d unmatched)", n_before, n_after, n_before - n_after))
  }
  message(sprintf("  %s distribution:", tgt))
  ct <- sort(table(data[[tgt]], useNA = "ifany"), decreasing = TRUE)
  for (i in seq_along(ct)) {
    message(sprintf("    %-30s %d", names(ct)[i], ct[i]))
  }

  data
}

################################################################################
# SECTION D: VALIDATION
################################################################################

#' Validate group sizes in a cross-tabulation
#'
#' @param data A dataframe
#' @param group_cols Character vector: column names to cross-tabulate
#' @param min_n Integer: minimum acceptable cell size
#' @return Invisible tibble with cross-tab and viability flags
validate_groups <- function(data, group_cols, min_n = 5) {
  stopifnot(all(group_cols %in% names(data)))

  counts <- data %>%
    count(across(all_of(group_cols)), name = "n") %>%
    arrange(across(all_of(group_cols)))

  counts$viable <- ifelse(counts$n >= min_n, "OK", "** LOW **")

  message(sprintf("\n=== Group Validation (min_n = %d) ===", min_n))
  message(sprintf("Total samples: %d", nrow(data)))

  # Print as wide table if exactly 2 grouping columns
  if (length(group_cols) == 2) {
    wide <- counts %>%
      mutate(display = sprintf("%d%s", n, ifelse(viable == "OK", "", " *"))) %>%
      select(-n, -viable) %>%
      pivot_wider(names_from = all_of(group_cols[2]),
                  values_from = display,
                  values_fill = "0 *")
    message("")
    print(as.data.frame(wide), row.names = FALSE)
  } else {
    message("")
    print(as.data.frame(counts), row.names = FALSE)
  }

  # Warnings for low cells
  low_cells <- counts %>% filter(viable == "** LOW **")
  if (nrow(low_cells) > 0) {
    message(sprintf("\nWARNING: %d group(s) below min_n = %d:", nrow(low_cells), min_n))
    for (i in seq_len(nrow(low_cells))) {
      row <- low_cells[i, ]
      label <- paste(group_cols, "=", row[group_cols], collapse = ", ")
      message(sprintf("  %s  (n=%d)", label, row$n))
    }
  } else {
    message(sprintf("\nAll groups meet minimum sample size (n >= %d).", min_n))
  }

  invisible(counts)
}

################################################################################
# SECTION E: DECISION SUMMARIES
################################################################################

#' Summarize what a category config excluded, collapsed, and kept
#'
#' Shows original counts for excluded/collapsed levels so you can verify
#' your decisions. Returns a list of tibbles suitable for kable() in Rmd.
#'
#' @param data_before Dataframe BEFORE applying the config (raw data)
#' @param config A category config list
#' @return Invisible list with tibbles: $kept, $collapsed, $excluded, $unaccounted
summarize_category_decisions <- function(data_before, config) {
  src <- config$source_col
  observed <- data_before %>%
    count(!!sym(src), name = "n") %>%
    arrange(desc(n)) %>%
    mutate(pct = round(n / sum(n) * 100, 1))

  accounted <- c(config$keep, unlist(config$collapse), config$exclude)
  unaccounted_levels <- setdiff(observed[[src]][!is.na(observed[[src]])], accounted)

  # Kept levels
  kept <- observed %>%
    filter(!!sym(src) %in% config$keep) %>%
    mutate(action = "Kept as-is", maps_to = !!sym(src))

  # Collapsed levels
  collapsed_rows <- list()
  for (target in names(config$collapse)) {
    for (source_level in config$collapse[[target]]) {
      row <- observed %>% filter(!!sym(src) == source_level)
      if (nrow(row) > 0) {
        collapsed_rows[[length(collapsed_rows) + 1]] <- row %>%
          mutate(action = "Collapsed", maps_to = target)
      }
    }
  }
  collapsed <- if (length(collapsed_rows) > 0) bind_rows(collapsed_rows) else
    tibble(!!sym(src) := character(), n = integer(), pct = double(),
           action = character(), maps_to = character())

  # Excluded levels
  excluded <- observed %>%
    filter(!!sym(src) %in% config$exclude) %>%
    mutate(action = "Excluded", maps_to = NA_character_)

  # Unaccounted levels
  unacc <- observed %>%
    filter(!!sym(src) %in% unaccounted_levels) %>%
    mutate(action = "Unaccounted (dropped)", maps_to = NA_character_)

  # Combined summary
  all_decisions <- bind_rows(kept, collapsed, excluded, unacc) %>%
    select(!!sym(src), n, pct, action, maps_to) %>%
    arrange(desc(n))

  result <- list(
    all = all_decisions,
    kept = kept,
    collapsed = collapsed,
    excluded = excluded,
    unaccounted = unacc,
    total_before = nrow(data_before),
    total_after = sum(kept$n) + sum(collapsed$n)
  )

  invisible(result)
}

#' Summarize what a grouping config matched and excluded
#'
#' Shows which rows matched each rule (with counts from source columns)
#' and what was excluded as "Other". Returns a list of tibbles.
#'
#' @param data_before Dataframe BEFORE applying the config
#' @param config A grouping config list
#' @return Invisible list with tibbles: $rule_matches, $unmatched
summarize_grouping_decisions <- function(data_before, config) {
  # Temporarily apply to get assignments without modifying original
  temp <- data_before
  tgt <- config$target_col
  temp[[tgt]] <- config$other_label

  for (i in rev(seq_along(config$rules))) {
    rule_name <- names(config$rules)[i]
    rule <- config$rules[[i]]
    mask <- rep(TRUE, nrow(temp))
    for (col_name in names(rule)) {
      col_match <- temp[[col_name]] == rule[[col_name]]
      col_match[is.na(col_match)] <- FALSE
      mask <- mask & col_match
    }
    temp[[tgt]][mask] <- rule_name
  }

  # Identify source columns used in rules
  rule_cols <- unique(unlist(lapply(config$rules, names)))

  # Matched rows summary
  matched <- temp %>%
    filter(!!sym(tgt) != config$other_label) %>%
    count(across(all_of(c(rule_cols, tgt))), name = "n") %>%
    arrange(!!sym(tgt), desc(n))

  # Unmatched rows summary
  unmatched <- temp %>%
    filter(!!sym(tgt) == config$other_label) %>%
    count(across(all_of(rule_cols)), name = "n") %>%
    arrange(desc(n))

  # Summary by target group
  group_summary <- temp %>%
    count(!!sym(tgt), name = "n") %>%
    arrange(desc(n)) %>%
    mutate(pct = round(n / sum(n) * 100, 1))

  result <- list(
    group_summary = group_summary,
    rule_matches = matched,
    unmatched = unmatched,
    total_before = nrow(data_before),
    total_matched = sum(temp[[tgt]] != config$other_label),
    total_unmatched = sum(temp[[tgt]] == config$other_label)
  )

  invisible(result)
}

message("covariate_explorer.R loaded successfully.")
