suppressWarnings(suppressMessages(library(argparse)))
suppressWarnings(suppressMessages(library(YRUtils)))
suppressWarnings(suppressMessages(library(vroom)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(clusterProfiler)))
suppressWarnings(suppressMessages(library(AnnotationDbi)))

parser <- ArgumentParser(description = "GSEA analysis.")
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
parser$add_argument("--metadata_tsv",
    type = "character", nargs = 1,
    action = "store", required = TRUE,
    help = "path to gene metadata TSV file (in which GO annotation file is given)"
)
parser$add_argument("--logfc_column",
    type = "character", nargs = 1,
    action = "store", required = TRUE,
    help = "logFC column name"
)
parser$add_argument("--gene_id_column",
    type = "character", nargs = 1,
    action = "store", required = TRUE,
    help = "gene IDs column name"
)
parser$add_argument("--gene_id_type",
    type = "character", default = "SYMBOL",
    nargs = 1, action = "store",
    choices = c("SYMBOL", "ENTREZID", "ENSEMBL"),
    help = "output directory"
)

args <- parser$parse_args()
args[["local_files"]] <- file.path(args[["input_dir"]], args[["names"]])
file_transfer(args[["files"]], args[["local_files"]], ops = c("soft", "hard", "copy"))

metadata_df <- vroom(args[["metadata_tsv"]], col_names = c("key", "value"))
metadata_vec <- setNames(metadata_df[["value"]], metadata_df[["key"]])
orgdb <- loadDb(metadata_vec["orgdb"])

for (local_file in args[["local_files"]]) {
    df <- vroom(local_file) %>%
        dplyr::select(all_of(c(args[["lofc_column"]], args[["gene_id_column"]]))) %>%
        set_colnames(c("logfc", "gene_id")) %>%
        na.omit() %>%
        distinct() %>%
        group_by(gene_id) %>%
        slice_max(abs(logfc)) %>%
        slice_sample(n = 1)
    df <- df[order(df[["logfc"]], decreasing = TRUE, na.last = TRUE), ]
    gene_list <- setNames(df[["logfc"]], df[["gene_id"]])

    ego <- gseGO(
        geneList = gene_list,
        ont = "ALL",
        OrgDb = orgdb,
        keyType = args[["gene_id_type"]],
        minGSSize = 10,
        maxGSSize = 1000,
        pvalueCutoff = 0.25,
        pAdjustMethod = "BH"
    )

    saveRDS(ego, file = file.path(
        args[["output_dir"]],
        gsub("\\.[a-zA-Z0-9]+$", ".GSEA_GO.ALL.rds", basename(local_file))
    ))
    vroom_write(arrange(ego@result, ONTOLOGY, desc(NES)),
        file = file.path(
            args[["output_dir"]],
            gsub("\\.[a-zA-Z0-9]+$", ".GSEA_GO.ALL.tsv", basename(local_file))
        )
    )
}
