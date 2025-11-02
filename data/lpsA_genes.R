suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
    require(memoise)
})

.load_lpsA_genes_ <- function() {
    load("raw/lpsA/lpsA.rda")
    lpsA$genes %>% tibble::rownames_to_column("Cell_ID") %>% tibble
}

source("cache_db.R")

.load_lpsA_genes <- memoise(.load_lpsA_genes_,cache=cache_db)

lpsA_genes=.load_lpsA_genes()

rm(cache_db)
rm(.load_lpsA_genes)
rm(.load_lpsA_genes_)
