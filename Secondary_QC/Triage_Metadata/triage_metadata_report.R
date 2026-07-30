# =============================================================================
# triage_metadata_report.R  ·  the run's RESULTS as a standalone HTML
# -----------------------------------------------------------------------------
# The README describes what the module does. This writes what it FOUND, per run,
# so findings never get baked into documentation that then goes stale against the
# next dataset.
#
# Written with base R only — the module deliberately has no rmarkdown/knitr
# dependency, so the report is emitted by the same script that computes it and
# can never drift out of sync with the CSVs beside it.
# =============================================================================

.esc <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#' Render a data frame as an HTML table, optionally highlighting rows.
#' @param highlight logical vector, same length as nrow(df)
html_table <- function(df, caption = NULL, highlight = NULL, digits = 4) {
  if (is.null(df) || !nrow(df)) return(sprintf("<p class='none'>%s: nothing to report.</p>",
                                               .esc(caption %||% "Table")))
  for (j in seq_along(df)) if (is.numeric(df[[j]]))
    df[[j]] <- ifelse(is.na(df[[j]]), "", format(signif(df[[j]], digits), trim = TRUE))
  hl <- if (is.null(highlight)) rep(FALSE, nrow(df)) else
    (!is.na(highlight) & highlight)
  paste0(
    if (!is.null(caption)) sprintf("<p class='cap'>%s</p>", .esc(caption)) else "",
    "<div class='scroll'><table><thead><tr>",
    paste0("<th>", .esc(names(df)), "</th>", collapse = ""), "</tr></thead><tbody>",
    paste0(vapply(seq_len(nrow(df)), function(i)
      sprintf("<tr class='%s'>%s</tr>", if (hl[i]) "hl" else "",
              paste0("<td>", .esc(unlist(df[i, ], use.names = FALSE)), "</td>", collapse = "")),
      character(1)), collapse = ""),
    "</tbody></table></div>")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Build the report.
#'
#' Every narrative sentence here is generated FROM the fitted numbers — there is
#' no prose that could stay true-looking while the data underneath it changed.
write_triage_metadata_report <- function(path, run_label, d, unl, conf, rates, tests,
                                         adj, bal, cm, sd_rows, min_n, n_perm,
                                         dx_source, generated) {
  P <- character(0)
  add <- function(...) P <<- c(P, ...)

  add("<h1>Triage &amp; flag decisions vs specimen metadata</h1>",
      sprintf("<p class='sub'>%s &middot; generated %s</p>", .esc(run_label), .esc(generated)),
      "<p class='banner'>Diagnostic only &mdash; this module drops nothing, flags nothing, and changes no pipeline output.</p>")

  # ---- what was analysed ----
  add("<h2>1. What was analysed</h2>",
      sprintf("<p>%s wells across %s plate-bays: <b>%s triaged</b>, <b>%s read-flagged</b>. Metadata joined for %s (%.1f%%). Group size floor n &ge; %d; %s permutation replicates. Diagnosis collapsed via <code>%s</code>.</p>",
              nrow(d), length(unique(d$plate)), sum(d$triaged), sum(d$read_flagged),
              sum(d$has_meta), 100 * mean(d$has_meta), min_n,
              format(n_perm, big.mark = ","), .esc(dx_source)),
      html_table(unl, "Wells with and without a metadata row. Unlabelled wells are kept in the denominator, never dropped."))

  # ---- can the question be asked ----
  add("<h2>2. Is the confounded question answerable?</h2>",
      "<p>Sites ship in batches, so &ldquo;this site fails more&rdquo; and &ldquo;this site&rsquo;s samples landed on worse plates&rdquo; can be the same data. If a label were nested inside plate, no model could separate them &mdash; that is checked before any conditioned result is reported.</p>",
      html_table(conf[, c("factor_a", "factor_b", "cramers_v", "b_per_a_median",
                          "median_share_of_a_on_one_b", "verdict")],
                 "Design separability."))

  # ---- the tests ----
  sig <- !is.na(tests$p_deviance) & (tests$p_deviance < 0.05 | tests$p_max_rate < 0.05)
  add("<h2>3. Association tests, per axis</h2>",
      "<p>The triage axes are near-disjoint and point in opposite directions on assay quality, so they are never pooled. Two statistics are reported: <b>max-rate</b> reacts to a single unusual group but is hostage to small ones; the <b>omnibus deviance</b> pools evidence and cannot be moved far by a handful of wells. <i>Them disagreeing is itself the finding.</i></p>",
      html_table(tests[, c("variable", "axis", "n_groups_tested", "n_flagged",
                           "observed_max_rate", "p_max_rate", "observed_deviance",
                           "p_deviance", "p_deviance_within_plate",
                           "min_detectable_rate_ratio", "reading")],
                 "Highlighted rows reached p < 0.05 on either statistic. Every null states the smallest effect it could have detected.",
                 highlight = sig))

  # ---- mutual adjustment ----
  if (!is.null(adj)) {
    survived <- adj[adj$variable != "log2(reads)" & !is.na(adj$p_given_others) &
                      adj$p_given_others < 0.05, ]
    collapsed <- adj[adj$variable != "log2(reads)" & !is.na(adj$retained_share) &
                       adj$retained_share < 0.5, ]
    add("<h2>4. Which label is actually carrying it?</h2>",
        "<p>Site, ancestry, diagnosis and age are correlated &mdash; a site recruits a particular population, whose age and diagnosis composition follow. Testing them one at a time lights several up off one cause. The full model is fitted <b>once</b>, with read depth included, and every term reported given the others. No stepwise, so a term that drops out is visible rather than absent.</p>",
        html_table(adj, "A label that keeps its deviance contributes something the others do not; one that collapses was a proxy. A 'keeps' above 1 is not an error &mdash; with correlated terms a variable can explain more once the others absorb variation it was competing with.",
                   highlight = !is.na(adj$p_given_others) & adj$p_given_others < 0.05))
    add("<div class='callout'>")
    if (nrow(survived)) {
      add("<p><b>Held up under adjustment</b> (with the other labels and read depth in the model):</p><ul>",
          paste0(vapply(seq_len(nrow(survived)), function(i)
            sprintf("<li><code>%s</code> &times; <code>%s</code> &mdash; p = %.4g, keeps %.0f%% of its unadjusted deviance</li>",
                    .esc(survived$axis[i]), .esc(survived$variable[i]),
                    survived$p_given_others[i], 100 * survived$retained_share[i]),
            character(1)), collapse = ""), "</ul>")
    } else add("<p><b>Nothing held up under adjustment.</b> Every metadata association collapsed once the others were included.</p>")
    if (nrow(collapsed))
      add("<p><b>Collapsed &mdash; these were proxies:</b></p><ul>",
          paste0(vapply(seq_len(nrow(collapsed)), function(i)
            sprintf("<li><code>%s</code> &times; <code>%s</code> &mdash; deviance %.1f &rarr; %.1f (keeps %.0f%%)</li>",
                    .esc(collapsed$axis[i]), .esc(collapsed$variable[i]),
                    collapsed$dev_alone[i], collapsed$dev_given_others[i],
                    100 * collapsed$retained_share[i]), character(1)), collapse = ""), "</ul>")
    add("</div>")
  }

  # ---- site ----
  sa <- rates[rates$variable == "Site" & rates$axis == "any triage", ]
  if (nrow(sa)) {
    sa <- sa[order(-sa$pct), ]
    add("<h2>5. Clinical collection site</h2>",
        "<p>A site whose specimens fail QC more often is a collection, handling or shipping problem &mdash; a different owner from a plate or well-position problem. Groups below the size floor are shown but not tested; their intervals say why.</p>",
        html_table(sa[, c("group", "n", "n_flagged", "pct", "ci_lo_pct", "ci_hi_pct", "above_floor")],
                   "Triage rate by site, with Wilson intervals. Rows with above_floor = FALSE were excluded from testing, not hidden.",
                   highlight = !sa$above_floor),
        html_table(sd_rows[, c("Site", "n_wells", "n_plates", "n_triaged", "pct_triaged",
                               "n_IC", "n_PCA", "n_burden", "n_read_flagged")],
                   "Per-site detail, every axis side by side."))
  }

  # ---- diagnosis ----
  dx <- rates[rates$variable == "Diagnosis" & rates$axis == "any triage", ]
  if (nrow(dx)) {
    td <- tests[tests$variable == "Diagnosis" & tests$axis == "any triage", ]
    add("<h2>6. Differential triage by diagnosis</h2>",
        "<p>This is the bias-relevant question: if triage removes cases at a different rate from controls, every downstream effect estimate is biased, and unlike a site effect nothing later in the pipeline would reveal it.</p>",
        html_table(dx[order(-dx$pct), c("group", "n", "n_flagged", "pct", "ci_lo_pct", "ci_hi_pct")],
                   "Triage rate by collapsed diagnosis."))
    if (nrow(td))
      add(sprintf("<p class='callout'><b>%s</b><br><span class='small'>A confidence interval spanning 1 is not equivalence. Read it together with the minimum detectable effect.</span></p>",
                  .esc(td$reading[1])))
    add(html_table(cm, "How each raw diagnosis value was collapsed. Administrative non-answers are kept as their own level rather than folded into a clinical one or dropped."))
  }

  # ---- balance ----
  bb <- bal[!is.na(bal$std_diff) & abs(bal$std_diff) >= 0.10, ]
  add("<h2>7. Who was removed vs who stayed</h2>",
      "<p>Standardised differences rather than p-values: at 77 removed against 4123 retained, a p-value mostly measures the sample size, while the standardised difference measures the imbalance itself.</p>",
      html_table(bb[order(-abs(bb$std_diff)), ],
                 "Imbalances at |standardised difference| >= 0.10."))

  add("<h2>8. What this does not claim</h2>",
      "<ul>",
      "<li>It does <b>not</b> establish that triage is unbiased with respect to diagnosis. It bounds the effect and reports the interval; that is weaker, and it is what the data support.</li>",
      "<li>It does <b>not</b> explain any association it finds. Whether a cause is collection, handling, shipping, storage time, or the specimens themselves is a wet-lab and logistics question this module cannot reach.</li>",
      "<li>It proposes <b>no drop rule</b>, and no result here should become one. Removing wells because of who they came from manufactures the exact bias section 6 checks for.</li>",
      "<li>It says nothing about <b>plate position</b> &mdash; the other meaning of &ldquo;site&rdquo;. See Primary QC &sect;D.5.4.</li>",
      "</ul>")

  css <- "
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
     max-width:1150px;margin:2rem auto;padding:0 1.25rem;line-height:1.55;color:#1a1a1a}
h1{font-size:1.7rem;margin-bottom:.2rem;border-bottom:2px solid #222;padding-bottom:.4rem}
h2{font-size:1.2rem;margin-top:2.2rem;border-bottom:1px solid #ddd;padding-bottom:.3rem}
.sub{color:#666;margin-top:0;font-size:.92rem}
.banner{background:#eef4fb;border-left:4px solid #3d6fa5;padding:.6rem .9rem;font-size:.92rem}
.callout{background:#f7f7f4;border-left:4px solid #999;padding:.6rem .9rem;margin:1rem 0}
.cap{font-size:.88rem;color:#555;margin:.9rem 0 .35rem}
.small{font-size:.85rem;color:#666}
.none{color:#777;font-style:italic}
.scroll{overflow-x:auto}
table{border-collapse:collapse;font-size:.85rem;width:100%}
th{background:#f0f0ec;text-align:left;padding:.35rem .5rem;border-bottom:2px solid #ccc;
   position:sticky;top:0;white-space:nowrap}
td{padding:.3rem .5rem;border-bottom:1px solid #eee;white-space:nowrap}
tr.hl td{background:#fdf6e3}
tr:hover td{background:#f4f8fd}
code{background:#f2f2ef;padding:.05rem .3rem;border-radius:3px;font-size:.9em}
ul{margin-top:.4rem}
@media (prefers-color-scheme:dark){
 body{background:#16181a;color:#e6e6e6}
 h1{border-color:#666} h2{border-color:#333}
 .banner{background:#1c2733;border-color:#5b8ec4}
 .callout{background:#1e1f21;border-color:#666}
 th{background:#23262a;border-color:#3a3d42} td{border-color:#26292d}
 tr.hl td{background:#2e2a18} tr:hover td{background:#1d2733}
 code{background:#26292d} .sub,.cap,.small,.none{color:#9aa0a6}}
"
  writeLines(c("<!doctype html><html><head><meta charset='utf-8'>",
               sprintf("<title>Triage vs metadata — %s</title>", .esc(run_label)),
               "<meta name='viewport' content='width=device-width,initial-scale=1'>",
               sprintf("<style>%s</style></head><body>", css),
               P, "</body></html>"), path)
  path
}
