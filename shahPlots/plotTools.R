x11 = function (...) grDevices::x11(...,type='cairo')

require(ggnewscale)

plot_genomic_links<-function(sid,ct,boi,YMAX) {

    # boi == bins of interest

    ob=get_links(sid,ct)

    #
    # tcn plotting data.frame
    #

    dt=ob$cnv %>%
        left_join(lpsA_chromInfo) %>%
        filter(bin.id %in% boi) %>%
        select(bin.id,bin.chrom,arm,band,med.tcn)

    #
    # Pad ends
    #

    header=dt[1,]
    header$bin.id=header$bin.id-1
    header$med.tcn=2

    footer=dt[nrow(dt),]
    footer$bin.id=footer$bin.id+1
    footer$med.tcn=2

    dt=bind_rows(list(header,dt,footer)) %>% mutate(R=row_number())


    #
    # Put in gaps for X-axis
    #

    SW=6
    ss=dt %>% group_by(arm) %>% summarize(max=max(R)) %>% slice(-nrow(.))
    spacer=bind_rows(rep(list(ss),SW)) %>% arrange(arm) %>% group_by(max) %>%
        mutate(med.tcn=case_when(row_number()==1 ~ 2,row_number()==SW ~ -1, T ~ NA)) %>%
        mutate(R=max+(row_number()/(SW+1))) %>% ungroup %>% select(-max)

    bcls=c("black","grey60")
    bandColors=dt %>% distinct(arm) %>% mutate(band.color=ifelse(row_number()%%2==0,"A","B"))
    if(bandColors %>% filter(arm=="12q") %>% pull(band.color)=="A") {
        names(bcls)=c("A","B")
    } else {
        names(bcls)=c("B","A")
    }

    dp=spacer %>% bind_rows(dt) %>% arrange(R) %>% mutate(R=row_number()) %>% left_join(bandColors)

    fixR=dp %>% filter(med.tcn==-1) %>% pull(R)
    dp[fixR,]=dp[fixR,] %>% mutate(med.tcn=2,band.color=ifelse(band.color=="A","B","A"))

    xrng=range(dp$R)

    p0=dp %>%
        ggplot(aes(R,med.tcn)) +
            theme_light(12) +
            coord_cartesian(expand=F) +
            ylab("Median Copy Number") +
            xlab(NULL) +
            ggtitle(paste(sid,ct))

    tcn=dp %>% select(bin.id,bin.coor=R,med.tcn)

    links=ob$links %>%
        filter(Type != "A->A[!12]") %>%
        left_join(tcn,by=c(BIN1="bin.id")) %>%
        rename(bin.coor.1=bin.coor,tcn.1=med.tcn) %>%
        left_join(tcn,by=c(BIN2="bin.id")) %>%
        rename(bin.coor.2=bin.coor,tcn.2=med.tcn) %>%
        filter(!is.na(bin.coor.1) & !is.na(bin.coor.2)) %>%
        mutate(band.color=NA)

    rr=links$bin.coor.1 < links$bin.coor.2

    t=links$bin.coor.1[rr]
    links$bin.coor.1[rr]=links$bin.coor.2[rr]
    links$bin.coor.2[rr]=t

    t=links$tcn.1[rr]
    links$tcn.1[rr]=links$tcn.2[rr]
    links$tcn.2[rr]=t

    xaxis=dp %>% group_by(arm) %>% summarize(mid=(min(R)+max(R))/2)

    linkTypes=c("12->12", "12->A", "A->B")
    linkCols=c(
        RColorBrewer::brewer.pal(5,"Purples")[5],
        RColorBrewer::brewer.pal(5,"Oranges")[5],
        RColorBrewer::brewer.pal(5,"Greens")[5]
    )

    names(linkCols)=linkTypes

    if(is.null(YMAX)) {
        YMAX=1.25*max(tcn$med.tcn,na.rm=T)
    }


    p1=p0 +
        geom_curve(
            aes(x=bin.coor.1,xend=bin.coor.2,y=tcn.1,yend=tcn.2,color=Type),
            dat=links,curvature=.225,alpha=.4
        ) +
        scale_y_continuous(limits=c(0,YMAX*1.25)) +
        scale_x_continuous(limits=c(xrng[1]-5,xrng[2]+5),breaks=xaxis$mid,labels=xaxis$arm) +
        scale_color_manual(values=linkCols,breaks=names(linkCols)) +
        new_scale_color() +
        geom_step(aes(color=band.color,group=1),data=dp,alpha=.85,lwd=1) +
        scale_color_manual(values=bcls,breaks=names(bcls)) +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.position="none"
        )

    hmga2_binid=lpsA_gene.index %>% filter(hgnc.symbol=="HMGA2") %>% pull(bin.id)
    hmga2_R=dp %>% filter(bin.id==hmga2_binid) %>% pull(R)
    cat("HMGA2 coor",hmga2_binid,hmga2_R,"\n")
    p1 + geom_vline(xintercept=hmga2_R,alpha=.3333) + annotate("text",x=hmga2_R,y=0.5,label=" HMGA2",size=3,hjust=0)

}

plot_12q_links<-function(sid,ct,YMAX=NA,curv=-.45) {

    ob=get_links(sid,ct)

    binCoords=lpsA_chromInfo %>%
        filter(arm=="12q") %>%
        mutate(bin.coor=(bin.start+bin.end)/2e6) %>%
        select(bin.id,bin.coor)

    tcn=left_join(binCoords,ob$cnv)

    XMAX=105

    if(is.na(YMAX)) {
        YMAX=1.25*max(tcn$med.tcn,na.rm=T)
    }

    links=ob$links %>%
        filter(CHR1=="chr12" & CHR2=="chr12") %>%
        left_join(tcn,by=c(BIN1="bin.id")) %>%
        rename(bin.coor.1=bin.coor,tcn.1=med.tcn) %>%
        left_join(tcn,by=c(BIN2="bin.id")) %>%
        rename(bin.coor.2=bin.coor,tcn.2=med.tcn) %>%
        filter(!is.na(bin.coor.1) & !is.na(bin.coor.2)) %>%
        select(bin.coor.1,bin.coor.2,n,tcn.1,tcn.2) %>%
        filter(bin.coor.1<XMAX & bin.coor.2<XMAX)

    pg=tcn %>% ggplot(aes(bin.coor,med.tcn)) +
        theme_light(14) +
        geom_step(linewidth=1) +
        scale_x_continuous(limits=c(55,XMAX)) +
        scale_y_continuous(limits=c(0,YMAX)) +
        xlab("Chromosome 12 [Mbase]") +
        ylab("Median Copy Number") +
        ggtitle(paste(sid,ct)) +
        geom_curve(
            aes(x=bin.coor.1,xend=bin.coor.2,y=tcn.1,yend=tcn.2,color=log2(n)),
            dat=links,curvature=curv,alpha=.4) +
        scale_color_gradientn(colors=RColorBrewer::brewer.pal(5,"Purples")[-1]) +
        coord_cartesian(expand=F) +
        guides(color="none")

    pg

}

