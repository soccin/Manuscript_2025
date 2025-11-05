source("load_shah.R")
source("tools.R")

x11 = function (...) grDevices::x11(...,type='cairo')

# SampleID=
# cloneType="EDC"
# cloneType="Adv.Tumor"

TCN.FILT=2
PE.FILT=3 #10
SIZE=3e6
require(tidyverse)
require(patchwork)

data(lpsA_sig)
data(lpsA_Y)
data(lpsA_chromInfo)

# Order arms, bands correctly
#
lpsA_chromInfo=lpsA_chromInfo %>%
    mutate(band=factor(band,levels=unique(band))) %>%
    mutate(arm=factor(arm,levels=unique(arm)))

binCoords=lpsA_chromInfo %>%
        mutate(bin.coor=(bin.start+bin.end)/2e6) %>%
        select(bin.id,bin.coor)

data(lpsA_gene.index)
genes=c("CDK4","MDM2","JUN","TCF21","HMGA2")
genesI=lpsA_gene.index %>%
    filter(hgnc.symbol %in% genes) %>%
    select(gene=hgnc.symbol,chrom,start,end,arm,bin.id)

mapSamps=map_vec(strsplit(basename(fs::dir_ls("data/raw/maps/v2023",regex="maps_.*rds")),"_"),2)

BLSAMPS=FALSE
if(!BLSAMPS) {
    sampleIDs=lpsA_Y %>%
    distinct(sampleID) %>%
    filter(!grepl("^BL",sampleID)) %>%
        pull
    sampleIDs=intersect(sampleIDs,paste0(mapSamps,"P"))
    cloneTypes=c("EDC","Adv.Tumor")
    cloneForBins="Adv.Tumor"
} else {
    sampleIDs=lpsA_Y %>%
    distinct(sampleID) %>%
    filter(grepl("^BL",sampleID)) %>%
        pull
    sampleIDs=intersect(sampleIDs,paste0(mapSamps,"B"))
    cloneTypes=c("Stroma","Adv.Lipoma")
    cloneForBins="Adv.Lipoma"
}

sid="WD0539P"
ct="Adv.Tumor"

if(file.exists("YMAX_Intra12.csv")) {
  ymax=read_csv("YMAX_Intra12.csv")
  ymax=transpose(ymax)
  names(ymax)=map_vec(ymax,"sid")
}

pp=list()
sampleIDs=sid
halt("DDDDDD")
for(sid in sampleIDs) {
    cat("sid =",sid,"\n")

    #
    # Get regions of interest (ROI)
    # use adv.tumor for this
    #

    boi=get_bins_of_interest(sid,size=SIZE,cloneForBins)

    #
    # Scale y-axis the same for all clones
    #
    #dd=load_shah_allClones(sid)
    #cells=dd$Y %>% filter(clone.type %in% cloneTypes) %>% pull(cellID)
    #max_tcn=dd$X %>% filter(cellID %in% cells) %>% pull(tcn) %>% max

    for(ct in cloneTypes) {
        print(c(sid,ct))
        pp[[cc(sid,ct)]]=plot_genomic_links(sid,ct,boi,YMAX=ymax[[sid]]$YMAX) + xlab("Genome Pos")
    }
}

pdf(file=cc("shahPlot","Genomic",ifelse(BLSAMPS,"BL_samps","ALL"),TCN.FILT,PE.FILT,round(SIZE/1e6),"V2b.pdf"),width=14,height=5)
for(ii in seq(len(pp)/2)) {
    print(pp[[ii*2-1]]+pp[[ii*2]])
}
dev.off()

plt=list()
for(ii in seq(len(pp)/2)) {
    plt[[ii]]=pp[[ii*2-1]]+pp[[ii*2]]
}
file=cc("shahPlot","Genomic",ifelse(BLSAMPS,"BL_samps","ALL"),TCN.FILT,PE.FILT,round(SIZE/1e6),"V2b.rds")
saveRDS(plt,file=file,compress=T)

file=cc("shahPlot","Genomic",ifelse(BLSAMPS,"BL_samps","ALL"),TCN.FILT,PE.FILT,round(SIZE/1e6),"Separate_V2b.rds")
saveRDS(pp,file=file,compress=T)

