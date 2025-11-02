suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
    require(memoise)
})

.load_lpsA_chromInfo_ <- function() {
    load("raw/lpsA/lpsA.rda")
    tibble(lpsA$chromInfo) %>% select(1:9)
}

source("cache_db.R")

.load_lpsA_chromInfo <- memoise(.load_lpsA_chromInfo_,cache=cache_db)

lpsA_chromInfo=.load_lpsA_chromInfo()

rm(cache_db)
rm(.load_lpsA_chromInfo)
rm(.load_lpsA_chromInfo_)
