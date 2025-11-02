suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
    require(memoise)
})

.load_lpsA_sig_ <- function() {
    load("raw/lpsA/lpsA.rda")
    digest::digest(lpsA)
}

source("cache_db.R")

.load_lpsA_sig <- memoise(.load_lpsA_sig_,cache=cache_db)

lpsA_sig=.load_lpsA_sig()

rm(cache_db)
rm(.load_lpsA_sig)
rm(.load_lpsA_sig_)
