cache_db=cachem::cache_disk("./cache")

.load_shah<-function(SampleID,CloneType) {

    data(lpsA_Y)
    Y=lpsA_Y %>% filter(sampleID==SampleID & clone.type==CloneType)
    if(nrow(Y)==0) {

        cat("\n\tNo Samples matching",SampleID,"\n")
        cat("\tValid samples:",lpsA_Y %>% distinct(sampleID) %>% pull,"\n")
        cat("\n")
        quit()

    }

    data(lpsA_X)
    X=lpsA_X %>% filter(cellID %in% Y$cellID)
    maps=readRDS(file.path("../Circos",cc("maps",Y$bioID[1],".rds"))) %>%
        filter(cellID %in% Y$cellID) %>%
        filter(BIN1!=BIN2) %>%
        mutate(BIN1=gsub("B","",BIN1)%>%as.numeric,BIN2=gsub("B","",BIN2)%>%as.numeric) %>%
        filter(abs(BIN1-BIN2)>1) %>%
        group_by(CHR1,BIN1,CHR2,BIN2,N,Nx,cellID) %>%
        summarize(n=sum(n)) %>%
        ungroup

    ret=list(X=X,Y=Y,maps=maps)

}

load_shah <- memoise::memoise(.load_shah,cache=cache_db)

.load_shah_allClones<-function(SampleID) {

    data(lpsA_Y)
    Y=lpsA_Y %>% filter(sampleID==SampleID & !is.na(ConsensusC))
    if(nrow(Y)==0) {

        cat("\n\tNo Samples matching",SampleID,"\n")
        cat("\tValid samples:",lpsA_Y %>% distinct(sampleID) %>% pull,"\n")
        cat("\n")
        quit()

    }

    data(lpsA_X)
    X=lpsA_X %>% filter(cellID %in% Y$cellID)
    maps=readRDS(file.path("../Circos",cc("maps",Y$bioID[1],".rds"))) %>%
        filter(cellID %in% Y$cellID) %>%
        filter(BIN1!=BIN2) %>%
        mutate(BIN1=gsub("B","",BIN1)%>%as.numeric,BIN2=gsub("B","",BIN2)%>%as.numeric) %>%
        group_by(CHR1,BIN1,CHR2,BIN2,N,Nx,cellID) %>%
        summarize(n=sum(n)) %>%
        ungroup %>%
        left_join(lpsA_Y %>% select(cellID,ConsensusC,clone.type))

    ret=list(X=X,Y=Y,maps=maps)

}

load_shah_allClones <- memoise::memoise(.load_shah_allClones,cache=cache_db)
