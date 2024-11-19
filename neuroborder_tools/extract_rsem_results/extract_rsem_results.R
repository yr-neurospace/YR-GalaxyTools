suppressWarnings(suppressMessages(library(argparse)))
suppressWarnings(suppressMessages(library(YRUtils)))
suppressWarnings(suppressMessages(library(vroom)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(magrittr)))

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

parser <- ArgumentParser(description = "Extract RSEM results.")
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

args <- parser$parse_args()
args[["local_files"]] <- file.path(args[["input_dir"]], args[["names"]])
file_transfer(args[["files"]], args[["local_files"]], ops = c("soft", "hard", "copy"))

if (!all(str_detect(basename(args[["local_files"]]), CONST_LS[[args[["type"]]]]["file_pattern"]))) {
    stop("invalid file name format")
}

# generate sample sheet
sample_df <- tibble(
    file = args[["local_files"]],
    basename = basename(file)
)
sample_df <- bind_cols(sample_df, basename(sample_df[["basename"]]) %>%
    lapply(function(x) {
        str_split(x, fixed("_"))[[1]][1:2]
    }) %>%
    do.call(rbind, .) %>% as.data.frame() %>%
    set_colnames(c("group", "replicate")) %>%
    mutate(sample = paste0(group, "_", replicate)))

vroom_write(sample_df,
    file = file.path(args[["output_dir"]], "sample_sheet.tsv"),
    delim = "\t", append = FALSE, col_names = TRUE
)

# extract raw counts, TPMs, and FPKMs
metadata_df <- vroom(args[["metadata_tsv"]], col_names = c("key", "value"))
metadata_vec <- setNames(metadata_df[["value"]], metadata_df[["key"]])
mpt <- vroom(metadata_vec[CONST_LS[[args[["type"]]]]["metadata_field"]])

expr_cols <- c("expected_count", "TPM", "FPKM")
for (expr_type in expr_cols) {
    df <- tibble()
    for (i in seq_len(nrow(sample_df))) {
        df <- vroom(sample_df[i, ][["file"]]) %>%
            select(all_of(c(setNames(CONST_LS[[args[["type"]]]]["id_field"], NULL), expr_type))) %>%
            mutate(sample = sample_df[i, ][["sample"]]) %>%
            distinct() %>%
            bind_rows(df)
    }
    df <- pivot_wider(df,
        id_cols = all_of(setNames(CONST_LS[[args[["type"]]]]["id_field"], NULL)),
        names_from = "sample", values_from = all_of(expr_type), values_fill = 0, names_sort = TRUE
    )
    df <- right_join(mpt, df, by = setNames(CONST_LS[[args[["type"]]]]["id_field"], NULL)) %>%
        distinct()
    vroom_write(df, file = file.path(args[["output_dir"]], paste0(args[["type"]], ".", expr_type, ".tsv")))
}
