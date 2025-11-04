get_bins_of_interest<-function(sid,size=100e6,ct) {

    ob=get_links(sid,ct)

    linkBins=ob$links %>% filter(CHR1!=CHR2) %>% select(BIN1,BIN2) %>% unlist %>% unname

    band.rois=ob$cnv %>%
        left_join(lpsA_chromInfo) %>%
        arrange(desc(med.tcn)) %>%
        filter(arm!="12q") %>%
        distinct(band,.keep_all=T) %>%
        mutate(linked=bin.id %in% linkBins) %>%
        filter(linked) %>%
        mutate(cum.width=cumsum(width)) %>%
        filter(cum.width<=size) %>%
        pull(band)

    lpsA_chromInfo %>%
        filter(band %in% band.rois | (arm=="12q" & bin.start>55e6 & bin.end<105e6) ) %>%
        pull(bin.id)

}

get_links<-function(sid,ct) {

    cat("sid =",sid,"ct =",ct,"\n")

    oo=load_shah(sid,ct)

    #
    # Normalize for cell count by scaling
    # to 100 fixed cells
    #
    cellNumScaling=lpsA_Y %>%
        filter(sampleID==sid & clone.type==ct) %>%
        count() %>%
        mutate(scaleCells=100/n) %>%
        pull

    maps=oo$maps %>%
        mutate(NC=cellNumScaling*n/(N/1e6)) %>%
        group_by(CHR1,BIN1,CHR2,BIN2) %>%
        summarize(n=sum(NC)) %>%
        ungroup %>%
        mutate(Type=case_when(
            CHR1=="chr12" & CHR2=="chr12" ~ "12->12",
            CHR1=="chr12" | CHR2=="chr12" ~ "12->A",
            CHR1!=CHR2 ~ "A->B",
            T ~ "A->A[!12]"
        ))

    cnv=oo$X %>%
        mutate(bin.id=as.numeric(bin.id)) %>%
        group_by(bin.id) %>%
        summarize(med.tcn=median(tcn))

    tcnFilterBins=cnv %>% filter(med.tcn>TCN.FILT) %>% pull(bin.id)

    links=maps %>%
        filter(n>PE.FILT) %>%
        filter(BIN1 %in% tcnFilterBins & BIN2 %in% tcnFilterBins)

    list(links=links,cnv=cnv)

}
