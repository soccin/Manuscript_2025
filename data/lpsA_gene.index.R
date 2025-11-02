suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
    require(memoise)
})

.load_lpsA_gene.index_ <- function() {
    load("raw/lpsA/lpsA.rda")
    lpsA$gene.index %>% tibble::rownames_to_column("Cell_ID") %>% tibble
}

source("cache_db.R")

.load_lpsA_gene.index <- memoise(.load_lpsA_gene.index_,cache=cache_db)

lpsA_gene.index=.load_lpsA_gene.index()

rm(cache_db)
rm(.load_lpsA_gene.index)
rm(.load_lpsA_gene.index_)
