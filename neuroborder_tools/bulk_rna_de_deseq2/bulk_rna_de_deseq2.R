suppressWarnings(suppressMessages(library(argparse)))
suppressWarnings(suppressMessages(library(YRUtils)))
suppressWarnings(suppressMessages(library(vroom)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(tximport)))
suppressWarnings(suppressMessages(library(DESeq2)))

CONST_LS <- list(
    gene = c(
        "file_pattern" = "^[a-zA-Z]+[a-zA-Z0-9]*_rep[0-9]+_anno_rsem\\.genes\\.results$",
        "id_field" = "gene_id",
        "metadata_field" = "gene_id_name_mapping_table"
    ),
    isoform = c(
        "file_pattern" = "^[a-zA-Z]+[a-zA-Z0-9]*_rep[0-9]+_anno_rsem\\.isoforms\\.results$",
        "id_field" = "transcript_id",
        "metadata_field" = "transcript_id_name_mapping_table"
    )
)

parser <- ArgumentParser(description = "Bulk RNA-seq differential expression analysis.")
parser$add_argument("--names",
    type = "character", nargs = "+",
    action = "extend", required = TRUE,
    help = "input file names separated by white space"
)
parser$add_argument("--files",
    type = "character", nargs = "+",
    action = "extend", required = TRUE,
    help = "input file paths separated by white space"
)
parser$add_argument("--input_dir",
    type = "character", nargs = 1,
    action = "store", required = TRUE,
    help = "input directory"
)
parser$add_argument("--output_dir",
    type = "character", nargs = 1,
    action = "store", required = TRUE,
    help = "output directory"
)
parser$add_argument("--type",
    type = "character", default = "gene",
    nargs = 1, action = "store",
    choices = c("gene", "isoform"),
    help = "output directory"
)
parser$add_argument("--metadata_tsv",
    type = "character", nargs = 1,
    action = "store", required = TRUE,
    help = "path to gene metadata TSV file (in which gene/transcript IDs to names mapping table are given)"
)
parser$add_argument("--sample_sheet_tsv",
    type = "character", nargs = 1,
    action = "store", required = TRUE,
    help = "path to sample sheet TSV file"
)
parser$add_argument("--tpm_tsv",
    type = "character", nargs = 1,
    action = "store", required = TRUE,
    help = "path to gene/isoform TPM TSV file (used to filter DEGs)"
)
parser$add_argument("--padj_th",
    type = "double", nargs = 1,
    action = "store", default = 0.05,
    help = "adjusted p-value threshold"
)
parser$add_argument("--logfc_th",
    type = "double", nargs = 1,
    action = "store", default = 1,
    help = "logFC threshold"
)
parser$add_argument("--tpm_th",
    type = "double", nargs = 1,
    action = "store", default = 1,
    help = "TPM/FPKM threshold"
)

args <- parser$parse_args()
args[["local_files"]] <- file.path(args[["input_dir"]], args[["names"]])
file_transfer(args[["files"]], args[["local_files"]], ops = c("soft", "hard", "copy"))

if (!all(str_detect(basename(args[["local_files"]]), CONST_LS[[args[["type"]]]]["file_pattern"]))) {
    stop("invalid file name format")
}

# read in sample sheet
sample_sheet_df <- vroom(args[["sample_sheet_tsv"]]) %>%
    select(-any_of("file")) %>%
    as.data.frame()
sample_df <- inner_join(sample_sheet_df,
    data.frame(
        file = args[["local_files"]],
        basename = args[["names"]]
    ),
    by = "basename"
) %>%
    arrange(group, replicate) %>%
    mutate(group = factor(group))
row.names(sample_df) <- sample_df[["sample"]]
if (nrow(sample_df) != nrow(sample_sheet_df)) {
    stop("input files are not exactly equivalent to those files within sample sheet file")
}

# read in quantification files
if (args[["type"]] == "gene") {
    tx.rsem <- tximport(setNames(sample_df[["file"]], sample_df[["sample"]]), type = "rsem", txIn = FALSE, txOut = FALSE)
} else if (args[["type"]] == "isoform") {
    tx.rsem <- tximport(setNames(sample_df[["file"]], sample_df[["sample"]]), type = "rsem", txIn = TRUE, txOut = TRUE)
}

# filter out genes/transcripts with lengths <= 0
non_zero_length <- apply(tx.rsem$length, 1, function(x) {
    all(x > 0)
})

tx.rsem$abundance <- tx.rsem$abundance[non_zero_length, ]
tx.rsem$counts <- tx.rsem$counts[non_zero_length, ]
tx.rsem$length <- tx.rsem$length[non_zero_length, ]

# create DESeqDatatSet
ddsTxi <- DESeqDataSetFromTximport(
    txi = tx.rsem,
    colData = sample_df,
    design = ~group
)

# keep rows that have at least N reads in total (this is just a simple filtering process and is optional)
keep <- rowSums(counts(ddsTxi)) >= 10
ddsTxi <- ddsTxi[keep, ]

group_levels <- levels(sample_df[["group"]])
de_ls <- list(raw = list(), na_filter = list(), unique_symbols = list(), only_degs = list())
for (ref in group_levels) {
    compara_levels <- c(ref, group_levels[-which(group_levels == ref)])

    message(paste0("\n\nLevels: ", paste0(compara_levels, collapse = ", ")))
    message("Baseline: ", ref, "\n")

    # set the reference sample
    ddsTxi$group <- factor(ddsTxi$group, levels = compara_levels)
    # DE analysis
    ddsTxi <- DESeq(ddsTxi)

    # get results table
    res_tabs <- resultsNames(ddsTxi)
    res_tabs <- res_tabs[res_tabs != "Intercept"]
    message("\n", paste0(res_tabs, collapse = ", "))

    # save pairwise comparisons
    for (pair in res_tabs) {
        res <- results(ddsTxi, name = pair)
        res <- as.data.frame(res)
        res[[CONST_LS[[args[["type"]]]]["id_field"]]] <- row.names(res)
        res[["pvalue"]][res[["pvalue"]] == 0] <- .Machine[["double.xmin"]]
        res[["padj"]][res[["padj"]] == 0] <- .Machine[["double.xmin"]]
        de_ls[["raw"]][[gsub("^group_", "", pair)]] <- as_tibble(res)
    }
}

# attach gene info
metadata_df <- vroom(args[["metadata_tsv"]], col_names = c("key", "value"))
metadata_vec <- setNames(metadata_df[["value"]], metadata_df[["key"]])
mpt <- vroom(metadata_vec[CONST_LS[[args[["type"]]]]["metadata_field"]])
for (pair in names(de_ls[["raw"]])) {
    de_ls[["raw"]][[pair]] <- inner_join(mpt, de_ls[["raw"]][[pair]],
        by = setNames(CONST_LS[[args[["type"]]]]["id_field"], NULL)
    )
}

# filter rows containing NAs
for (pair in names(de_ls[["raw"]])) {
    pairs <- str_split(pair, fixed("_vs_"))[[1]]
    de_ls[["na_filter"]][[pair]] <- de_ls[["raw"]][[pair]] %>%
        filter(!(is.na(baseMean) | is.na(log2FoldChange) | is.na(pvalue) | is.na(padj))) %>%
        mutate(diff_flag = if_else(padj < args[["padj_th"]],
            if_else(abs(log2FoldChange) > args[["logfc_th"]],
                if_else(log2FoldChange > args[["logfc_th"]],
                    paste0(pairs[1], " Up"),
                    paste0(pairs[2], " Up")
                ),
                "NO"
            ),
            "NO"
        ))
}

if (args[["type"]] == "gene") {
    # filter duplicated gene symbols
    tpm <- vroom(args[["tpm_tsv"]])
    for (pair in names(de_ls[["na_filter"]])) {
        # filtered by TPMs
        pairs <- str_split(pair, fixed("_vs_"))[[1]]
        flag <- rep(FALSE, nrow(tpm))
        for (s in pairs) {
            tmp_df <- tpm[, str_detect(names(tpm), paste0("^", s, "_rep[0-9]+$"))]
            if (ncol(tmp_df) == 0) {
                stop("sample columns in TPM are not matched with those used to perform DE analysis")
            }
            flag <- flag | (rowSums(tmp_df > args[["tpm_th"]]) == ncol(tmp_df))
        }
        tmp_df <- de_ls[["na_filter"]][[pair]] %>%
            filter(gene_id %in% tpm[["gene_id"]][flag]) %>%
            mutate(gene_name = if_else(is.na(gene_name), gene_id, gene_name))

        # filtered by gene versions
        if ("gene_version" %in% names(tmp_df)) {
            tmp_df <- tmp_df %>%
                group_by(gene_name) %>%
                slice_max(gene_version) %>%
                ungroup()
        }

        # filtered by baseMean, padj, log2FoldChange
        # if still duplicated, sample one randomly
        de_ls[["unique_symbols"]][[pair]] <- tmp_df %>%
            group_by(gene_name) %>%
            slice_max(baseMean) %>%
            slice_min(padj) %>%
            slice_max(log2FoldChange) %>%
            slice_sample(n = 1) %>%
            ungroup()
    }
}

if (args[["type"]] == "isoform") {
    # filter duplicated isoform IDs
    tpm <- vroom(args[["tpm_tsv"]])
    for (pair in names(de_ls[["na_filter"]])) {
        # filtered by TPMs
        pairs <- str_split(pair, fixed("_vs_"))[[1]]
        flag <- rep(FALSE, nrow(tpm))
        for (s in pairs) {
            tmp_df <- tpm[, str_detect(names(tpm), paste0("^", s, "_rep[0-9]+$"))]
            if (ncol(tmp_df) == 0) {
                stop("sample columns in TPM are not matched with those used to perform DE analysis")
            }
            flag <- flag | (rowSums(tmp_df > args[["tpm_th"]]) == ncol(tmp_df))
        }
        tmp_df <- de_ls[["na_filter"]][[pair]] %>%
            filter(transcript_id %in% tpm[["transcript_id"]][flag])

        # filtered by transcript versions
        if ("transcript_version" %in% names(tmp_df)) {
            tmp_df <- tmp_df %>%
                group_by(transcript_id) %>%
                slice_max(transcript_version) %>%
                ungroup()
        }

        de_ls[["unique_symbols"]][[pair]] <- tmp_df
    }
}

# keep only DEGs
for (pair in names(de_ls[["unique_symbols"]])) {
    de_ls[["only_degs"]][[pair]] <- de_ls[["unique_symbols"]][[pair]] %>%
        filter(diff_flag != "NO")
}

# save results to files
setwd(args[["output_dir"]])
for (category in names(de_ls)) {
    dir.create(category)
    for (pair in names(de_ls[[category]])) {
        vroom_write(de_ls[[category]][[pair]], file = file.path(category, paste0(pair, ".tsv")))
    }

    cmd <- paste0("tar --use-compress-program=pigz -cvf ", category, ".tgz ", category)
    system(cmd, wait = TRUE)
}
