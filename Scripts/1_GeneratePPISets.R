#!/usr/bin/env Rscript
# ==========================================================
# NIRRA — Network-Informed Restricted-set Ridge Analysis
#
# Copyright 2026 William Tower
#
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
# ==========================================================
# ==========================================================
# DensePPI_v2.8 — All-hubs (FirstProtein→LastProtein) STRING subnetwork set generation (robust hub sampling + explicit evidence sourcing)
# ==========================================================
# Will Tower — 2026-03-18
# ==========================================================

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(igraph)
  library(tibble)
  library(readr)
})

# ==========================================================
# 1. Paths and parameters
# ==========================================================
in_file <- "~/path/to/inputmatrix.csv"
out_dir <- "~/path/to/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

string_api <- "https://string-db.org/api"
output_format <- "json"
method <- "network"
species <- 9606
confidence_cutoff <- 600

collection_min <- 2
collection_max <- 6
replicates_per_hub <- 10
ppi_alpha <- 0.05
autosave_every <- 500

cache_file <- file.path(out_dir, "string_cache_noTextMining.rds")
log_file <- file.path(out_dir, "DensePPI_missing_hubs.log")

set.seed(1)

# ==========================================================
# 2. Load dataset and define protein hubs
# ==========================================================
df <- read_csv(in_file, show_col_types = FALSE)

prot_start <- which(names(df) == "FirstProtein")
prot_end   <- which(names(df) == "LastProtein")

if (length(prot_start) == 0 | length(prot_end) == 0)
  stop("❌ Could not locate FirstProtein or LastProtein.")

prot_cols <- names(df)[prot_start:prot_end]

# ---- Clean names for STRING ----
prot_clean <- gsub("\\.\\.\\..*$", "", prot_cols)
prot_clean <- gsub("\\.\\..*$", "", prot_clean)
prot_clean <- trimws(prot_clean)

message("Proteome size: ", length(prot_cols))
message("Unique cleaned: ", length(unique(prot_clean)))

protein_map <- tibble(original = prot_cols, clean = prot_clean)

# ==========================================================
# 3. STRING cache (FIXED: no environment)
# ==========================================================
if (file.exists(cache_file)) {
  string_cache <- readRDS(cache_file)
} else {
  string_cache <- list()
}

get_neighbors <- function(protein) {
  
  if (!is.null(string_cache[[protein]]))
    return(string_cache[[protein]])
  
  url <- paste0(
    string_api, "/", output_format, "/", method,
    "?identifiers=", protein,
    "&species=", species,
    "&required_score=", confidence_cutoff
  )
  
  res <- try(httr::GET(url), silent = TRUE)
  if (inherits(res, "try-error") || res$status_code != 200) return(NULL)
  
  txt <- httr::content(res, "text", encoding = "UTF-8")
  if (is.null(txt) || txt == "" || grepl("Error", txt, ignore.case = TRUE))
    return(NULL)
  
  df <- tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0)
    return(NULL)
  
  # ==========================================================
  # FIX: Remove ONLY textmining-only edges
  # Keep edges if ANY other evidence exists
  # ==========================================================
  evidence_cols <- c("neighborhood", "fusion", "cooccurence",
                     "coexpression", "experimental", "database")
  
  if (all(c("textmining", evidence_cols) %in% names(df))) {
    other_evidence <- rowSums(df[, evidence_cols], na.rm = TRUE)
    df <- df[!(df$textmining > 0 & other_evidence == 0), , drop = FALSE]
  }
  
  if (nrow(df) == 0) return(NULL)
  
  tab <- tibble(
    protein1 = df$preferredName_A,
    protein2 = df$preferredName_B,
    combined_score = df$score * 1000
  )
  
  string_cache[[protein]] <<- tab
  saveRDS(string_cache, cache_file)
  
  return(tab)
}

# ==========================================================
# 4. Build STRING graph
# ==========================================================
message("Building STRING graph...")

unique_prots <- unique(prot_clean)
batch_size <- 500
batches <- split(unique_prots,
                 ceiling(seq_along(unique_prots)/batch_size))

links_list <- list()

for (b in seq_along(batches)) {
  
  batch <- batches[[b]]
  
  message("→ Batch ", b, "/", length(batches),
          " (", length(batch), ")")
  
  batch_links <- lapply(batch, get_neighbors)
  links_list[[b]] <- bind_rows(batch_links)
  
  Sys.sleep(1)  # safer pacing
}

links <- bind_rows(links_list) %>%
  filter(protein1 %in% prot_clean,
         protein2 %in% prot_clean) %>%
  distinct()

g <- igraph::graph_from_data_frame(links, directed = FALSE)

background_genes <- igraph::V(g)$name

message("Graph: ", length(background_genes), " nodes, ",
        igraph::gsize(g), " edges")

# ==========================================================
# 5. Hypergeometric enrichment (stable math)
# ==========================================================
get_ppi_p_hypergeom <- function(gene_set, graph, all_genes) {
  
  gs <- intersect(gene_set, all_genes)
  k <- length(gs)
  if (k < 2) return(NA_real_)
  
  subg <- igraph::induced_subgraph(graph, vids = gs)
  
  if (igraph::components(subg)$no > 1)
    return(NA_real_)
  
  x <- igraph::gsize(subg)
  
  N <- length(all_genes)
  M <- N * (N - 1) / 2
  
  K <- igraph::gsize(graph)
  m <- k * (k - 1) / 2
  
  stats::phyper(q = x - 1, m = K, n = M - K, k = m, lower.tail = FALSE)
}

# ==========================================================
# 6. Dense collection generation (FIXED HUB LOGIC)
# ==========================================================
sig_sets <- list()
seen_sets <- new.env(hash = TRUE, parent = emptyenv())

if (file.exists(log_file)) file.remove(log_file)

message("Generating dense sets...")

last_autosave <- 0

for (hub in background_genes) {
  
  neighbors <- igraph::V(g)[igraph::neighbors(g, hub)]$name
  
  # ==========================================================
  # FIX: allow degree-1 hubs
  # ==========================================================
  if (length(neighbors) == 0) {
    write(paste0("No edges for hub: ", hub),
          file = log_file, append = TRUE)
    next
  }
  
  found_any <- FALSE
  
  for (rep in seq_len(replicates_per_hub)) {
    
    for (k in collection_min:collection_max) {
      
      if (length(neighbors) < (k - 1)) next
      
      chosen <- unique(c(hub, sample(neighbors, k - 1)))
      
      subg <- igraph::induced_subgraph(g, vids = chosen)
      if (igraph::components(subg)$no > 1) next
      
      ppi_p <- get_ppi_p_hypergeom(chosen, g, background_genes)
      
      if (is.na(ppi_p) || ppi_p >= ppi_alpha) next
      
      found_any <- TRUE
      
      key <- paste(sort(chosen), collapse = ",")
      
      if (!exists(key, envir = seen_sets, inherits = FALSE)) {
        
        assign(key, TRUE, envir = seen_sets)
        
        sig_sets[[length(sig_sets) + 1]] <- tibble(
          Set       = length(sig_sets) + 1,
          Hub       = hub,
          Replicate = rep,
          Size      = k,
          PPI_p     = ppi_p,
          Proteins  = key
        )
      }
    }
  }
  
  if (!found_any) {
    write(paste0("No dense sets for hub: ", hub),
          file = log_file, append = TRUE)
  }
  
  # ==========================================================
  # FIX: clean autosave trigger
  # ==========================================================
  if ((length(sig_sets) - last_autosave) >= autosave_every) {
    
    message("💾 Autosaving at ", length(sig_sets), " sets")
    
    tmp <- bind_rows(sig_sets)
    
    write_csv(tmp,
              file.path(out_dir, "PPI_HighDensity_Collections_partial.csv")
    )
    
    last_autosave <- length(sig_sets)
  }
}

# ==========================================================
# 7. Final output
# ==========================================================
if (length(sig_sets) == 0) {
  message("⚠️ No significant sets generated.")
  quit(status = 0)
}

ppi_df <- bind_rows(sig_sets) %>%
  mutate(
    q_value = p.adjust(PPI_p, method = "fdr"),
    NegLog10PPI = -log10(PPI_p)
  ) %>%
  arrange(PPI_p)

write_csv(ppi_df,
          file.path(out_dir, "PPI_HighDensity_Collections_AllHubs.csv")
)

message("✅ Done: ", nrow(ppi_df), " unique dense sets")
message("📄 Log: ", log_file)