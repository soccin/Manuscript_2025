suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
    require(memoise)
})

.load_lpsA_XYZ_ <- function() {
    load("raw/lpsA/lpsA.rda")
    #
    # Transform data object you want
    # e.g.: lpsA$genes %>% tibble::rownames_to_column("Cell_ID") %>% tibble
    #
}

source("cache_db.R")

.load_lpsA_XYZ <- memoise(.load_lpsA_XYZ_,cache=cache_db)

lpsA_XYZ=.load_lpsA_XYZ()

rm(cache_db)
rm(.load_lpsA_XYZ)
rm(.load_lpsA_XYZ_)
