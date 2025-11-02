suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
    require(memoise)
})

.load_lpsA_X_ <- function() {
    load("raw/lpsA/lpsA.rda")
    tibble(lpsA$X) %>%
        tibble::rownames_to_column("bin.id") %>%
        gather(cellID,tcn,-bin.id)
}

source("cache_db.R")

.load_lpsA_X <- memoise(.load_lpsA_X_,cache=cache_db)

lpsA_X=.load_lpsA_X()

rm(cache_db)
rm(.load_lpsA_X)
rm(.load_lpsA_X_)
