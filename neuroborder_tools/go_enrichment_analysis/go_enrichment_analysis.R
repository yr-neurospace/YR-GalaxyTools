suppressWarnings(suppressMessages(library(argparse)))
suppressWarnings(suppressMessages(library(YRUtils)))
suppressWarnings(suppressMessages(library(vroom)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(clusterProfiler)))
suppressWarnings(suppressMessages(library(AnnotationDbi)))

parser <- ArgumentParser(description = "GO enrichment analysis.")
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
parser$add_argument("--group_column",
    type = "character", nargs = 1,
    action = "store", required = TRUE,
    help = "group column name"
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
        dplyr::select(all_of(c(args[["group_column"]], args[["gene_id_column"]]))) %>%
        set_colnames(c("group", "gene_id")) %>%
        na.omit() %>%
        distinct()
    gene_ls <- lapply(split(df, df[["group"]]), function(df) {
        unique(na.omit(df[["gene_id"]]))
    })

    if (length(gene_ls) > 0) {
        ego <- compareCluster(
            gene_ls,
            OrgDb = orgdb,
            keyType = args[["gene_id_type"]],
            fun = "enrichGO",
            ont = "ALL",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            minGSSize = 10,
            maxGSSize = 1000,
            readable = FALSE,
            pool = FALSE
        )

        if (!is.null(ego)) {
            saveRDS(ego, file = file.path(
                args[["output_dir"]],
                gsub("\\.[a-zA-Z0-9]+$", ".GO.ALL.rds", basename(local_file))
            ))
            vroom_write(arrange(ego@compareClusterResult, ONTOLOGY, Cluster),
                file = file.path(
                    args[["output_dir"]],
                    gsub("\\.[a-zA-Z0-9]+$", ".GO.ALL.tsv", basename(local_file))
                )
            )
        }
    }
}
