suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
    require(memoise)
})

.load_lpsA_Y_ <- function() {
    load("raw/lpsA/lpsA.rda")
    tibble(lpsA$Y)
}

source("cache_db.R")

.load_lpsA_Y <- memoise(.load_lpsA_Y_,cache=cache_db)

lpsA_Y=.load_lpsA_Y()

rm(cache_db)
rm(.load_lpsA_Y)
rm(.load_lpsA_Y_)
