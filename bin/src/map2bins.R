argv=commandArgs(trailing=T)

require(tidygenomics)
require(tidyverse)

bins=read_tsv(argv[1],col_names=F) %>% mutate(X2=X2+1) %>% mutate(BIN=sprintf("B%05d",row_number()))
map=read_tsv(argv[2],col_names=F)

if(nrow(map)==0) {
    cat("\nNo PE READS passed for",argv[2],"\n\n")
    quit()
}

map=map %>% mutate(X3=X2,X6=X5) %>% mutate(MID=sprintf("M%05d",row_number()))

bin1=genome_intersect(map,bins,by=c(X1="X1",X2="X2",X3="X3")) %>% select(MID,BIN1=BIN)
bin2=genome_intersect(map,bins,by=c(X4="X1",X5="X2",X6="X3")) %>% select(MID,BIN2=BIN)

map=map %>%
    left_join(bin1) %>%
    left_join(bin2) %>%
    filter(!is.na(BIN1) & !is.na(BIN2))

write_tsv(map,gsub(".bed.gz",".map.gz",argv[2]))

