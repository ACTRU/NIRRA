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
# ClusterSets_v2.2 — k-NN Jaccard Graph + Leiden Clustering for Large Collections (Correct deduplication)
# ==========================================================
# Will Tower — 2026-03-18
# ==========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(Matrix)
  library(data.table)
  library(igraph)
})

# ----------------------------
# 1) File paths
# ----------------------------
in_file <- "~/path/to/ALL_Set_Models.csv"
out_dir <- "~/path/to/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_csv        <- file.path(out_dir, "Full_Collections_Clustered_kNN_Jaccard.csv")
out_unique_csv <- file.path(out_dir, "Unique_Collections_Clustered_kNN_Jaccard.csv")
out_edges      <- file.path(out_dir, "OverlapGraph_kNN_Jaccard_Edges.csv")
out_metrics    <- file.path(out_dir, "OverlapGraph_kNN_Jaccard_RunMetrics.txt")

# ----------------------------
# 2) Parameters
# ----------------------------
k_neighbors  <- 200
min_overlap  <- 1
resolution   <- 0.01
use_leiden   <- TRUE

max_edges_pre_knn  <- Inf
max_edges_post_knn <- Inf

# ----------------------------
# 3) Helpers
# ----------------------------
stop_if_missing_col <- function(df, col) {
  if (!col %in% names(df)) stop("Missing required column: ", col)
}

clean_token <- function(x) {
  x <- gsub("\\.\\.\\.[0-9]+$", "", x)
  x <- gsub("\\.[0-9]+$", "", x)
  x
}

# 🔴 FIXED: canonicalize by sorting tokens
normalize_protein_string <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\\s+", "")
  
  toks <- strsplit(x, ",", fixed = TRUE)[[1]]
  toks <- toks[nchar(toks) > 0]
  
  toks <- clean_token(toks)
  toks <- unique(toks)
  toks <- sort(toks)  # ✅ critical fix
  
  paste(toks, collapse = ",")
}

write_metrics <- function(lines, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, con = path)
}

# ----------------------------
# 4) Load data
# ----------------------------
df <- readr::read_csv(in_file, show_col_types = FALSE)
message("✅ Loaded ", nrow(df), " total rows")

stop_if_missing_col(df, "Proteins")
df <- df %>% mutate(Proteins = as.character(Proteins))

n_rows_original <- nrow(df)

# ----------------------------
# 5) Normalize + collapse duplicates
# ----------------------------
message("🧬 Normalizing Proteins strings and collapsing duplicates...")

df$Proteins_Normalized <- vapply(df$Proteins, normalize_protein_string, character(1))

if (any(df$Proteins_Normalized == "" | is.na(df$Proteins_Normalized))) {
  bad <- which(df$Proteins_Normalized == "" | is.na(df$Proteins_Normalized))
  stop("Found ", length(bad), " invalid normalized Proteins. Example rows: ",
       paste(head(bad, 20), collapse = ", "),
       if (length(bad) > 20) " ..." else "")
}

dup_map <- data.frame(
  OriginalRow = seq_len(nrow(df)),
  Proteins_Normalized = df$Proteins_Normalized,
  stringsAsFactors = FALSE
)

df_unique <- df %>%
  distinct(Proteins_Normalized, .keep_all = TRUE)

n_unique <- nrow(df_unique)
n_dup_removed <- n_rows_original - n_unique

message("✅ Unique Protein rows: ", n_unique)
message("♻️ Duplicate rows collapsed: ", n_dup_removed)

df_unique$Proteins <- df_unique$Proteins_Normalized
n_sets <- nrow(df_unique)

# ----------------------------
# 6) Parse tokens
# ----------------------------
message("🧬 Parsing unique Proteins tokens...")

proteins_list <- strsplit(df_unique$Proteins, ",", fixed = TRUE)

proteins_list <- lapply(proteins_list, function(v) {
  v <- v[nchar(v) > 0]
  unique(v)
})

if (any(lengths(proteins_list) == 0)) {
  bad <- which(lengths(proteins_list) == 0)
  stop("Found empty Protein lists after parsing. Rows: ",
       paste(head(bad, 20), collapse = ", "))
}

all_proteins <- unique(unlist(proteins_list, use.names = FALSE))
message("🧬 Unique proteins: ", length(all_proteins))

# ----------------------------
# 7) Sparse incidence matrix
# ----------------------------
message("🧮 Building sparse incidence matrix...")

lens <- lengths(proteins_list)
i_idx <- rep.int(seq_along(proteins_list), times = lens)
j_idx <- match(unlist(proteins_list, use.names = FALSE), all_proteins)

mat <- Matrix::sparseMatrix(
  i = i_idx,
  j = j_idx,
  x = 1L,
  dims = c(n_sets, length(all_proteins)),
  giveCsparse = TRUE
)

set_sizes <- Matrix::rowSums(mat)

# ----------------------------
# 8) Overlap matrix
# ----------------------------
message("⚙️ Computing overlap matrix...")

overlap <- Matrix::tcrossprod(mat)
Matrix::diag(overlap) <- 0

overlap@x[overlap@x < min_overlap] <- 0
overlap <- Matrix::drop0(overlap)

if (length(overlap@x) == 0) {
  stop("No overlaps found — check token normalization.")
}

# ----------------------------
# 9) Extract edges
# ----------------------------
edges <- data.table::as.data.table(Matrix::summary(overlap))
setnames(edges, c("i", "j", "intersect"))

if (is.finite(max_edges_pre_knn) && nrow(edges) > max_edges_pre_knn) {
  setorder(edges, -intersect)
  edges <- edges[1:max_edges_pre_knn]
}

# ----------------------------
# 10) Jaccard + kNN
# ----------------------------
message("🧠 Computing Jaccard kNN...")

edges[, jaccard := intersect / (set_sizes[i] + set_sizes[j] - intersect)]

setorder(edges, i, -jaccard, -intersect, j)
edges_knn_i <- edges[, head(.SD, k_neighbors), by = i]

setorder(edges, j, -jaccard, -intersect, i)
edges_knn_j <- edges[, head(.SD, k_neighbors), by = j]

edges_knn <- unique(rbind(
  edges_knn_i[, .(i, j, intersect, jaccard)],
  edges_knn_j[, .(i, j, intersect, jaccard)]
))

rm(edges, overlap)
gc()

# ----------------------------
# 11) Graph + clustering
# ----------------------------
message("🧩 Building graph...")

g <- graph_from_data_frame(
  data.frame(
    from = as.character(edges_knn$i),
    to   = as.character(edges_knn$j),
    weight = edges_knn$jaccard
  ),
  directed = FALSE,
  vertices = data.frame(name = as.character(seq_len(n_sets)))
)

g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

if (use_leiden && "cluster_leiden" %in% getNamespaceExports("igraph")) {
  cl <- cluster_leiden(g, weights = E(g)$weight, resolution_parameter = resolution)
} else {
  cl <- cluster_louvain(g, weights = E(g)$weight)
}

df_unique$Cluster <- as.integer(membership(cl))

# ----------------------------
# 12) Map back to original rows
# ----------------------------
cluster_lookup <- df_unique %>%
  select(Proteins_Normalized, Cluster)

df_out <- df %>%
  left_join(cluster_lookup, by = "Proteins_Normalized")

# ----------------------------
# 13) Save outputs
# ----------------------------
write_csv(df_unique, out_unique_csv)
write_csv(df_out, out_csv)

# ----------------------------
# 14) Metrics
# ----------------------------
write_metrics(c(
  paste0("Original rows: ", n_rows_original),
  paste0("Unique rows: ", n_unique),
  paste0("Duplicates removed: ", n_dup_removed)
), out_metrics)

message("✔ COMPLETE — canonicalized clustering complete")
