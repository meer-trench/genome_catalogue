library(tidyverse)
library(dplyr)
library(networkD3)


setwd('/Users/jiahuizhu/Git/meer_sediments_metagenomics/Figure_Script/MEERv2.1/Taxonomy')
summary <- read_tsv("./MEERv21_gtdbr220_taxonomy.tsv")
summary %>% filter(Phylum=="Actinomycetota") %>% count()
unclassified_species<-summary %>% filter(is.na(Species)) %>% select(Kingdom,Phylum,Genus,Species)
alluvial_unknown <- unclassified_species %>% mutate(status=case_when(
  is.na(Genus) ~ "Genus+",
  .default = "Species",
)) %>% select(Kingdom,Phylum,status)
alluvial_unknown[which(is.na(alluvial_unknown$Phylum)),"Phylum"] <- "Unclassified"
alluvial_freq <- alluvial_unknown%>% group_by(Kingdom,Phylum,status) %>%  count()
bacteria_alluvial <- alluvial_unknown %>% group_by(Kingdom,Phylum) %>%  count() %>%filter(Kingdom=="Bacteria")%>%arrange(desc(n))
bac_top15 <- bacteria_alluvial[1:15,"Phylum"]
archaea_alluvial <- alluvial_unknown %>% group_by(Kingdom,Phylum) %>%  count() %>%filter(Kingdom=="Archaea")%>%arrange(desc(n))
arc_top <- archaea_alluvial[1:7,"Phylum"]
phylum_level <- arc_top %>% bind_rows(bac_top15)
plot_allu_feq <- alluvial_freq %>% filter(Phylum %in% as.matrix(phylum_level))

source_target1 <- alluvial_unknown %>% group_by(Kingdom,Phylum) %>%  count() %>% filter(Phylum %in% as.matrix(phylum_level))%>% as.data.frame()
colnames(source_target1) <- c("source","target","value")
source_target1 <- source_target1 %>% arrange(desc(value))
source_target2 <- alluvial_unknown %>% group_by(Phylum,status) %>%  count() %>% filter(Phylum %in% as.matrix(phylum_level))%>%as.data.frame()
colnames(source_target2) <- c("source","target","value")
source_target2 <- source_target2 %>% arrange(desc(value),desc(target))

source_target <- source_target1 %>% bind_rows(source_target2) %>% as.data.frame()
colnames(source_target) <- c("source","target","value")

nodes <- data.frame(
  name=c(as.character(source_target$source), 
         as.character(source_target$target)) %>% unique()
)
source_target$IDsource <- match(source_target$source, nodes$name)-1 
source_target$IDtarget <- match(source_target$target, nodes$name)-1
nw <- sankeyNetwork(Links = source_target, Nodes = nodes,
                    Source = "IDsource", Target = "IDtarget",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE,
                    width = 500,
                    height = 500,
                    fontSize=12,nodeWidth=10,nodePadding = 6,
                    fontFamily="Arial")
saveNetwork(nw,"./Prokayotes_Classification_genus+_Sankey.html")

unclassified_species$status <-""
unclassified_species$status[which(unclassified_species$Genus!="")] <- "Genus"
unclassified_species$status[which(unclassified_species$Genus=="")] <- "Family+"
alluvial_unknown <- unclassified_species[,c("Kingdom","Phylum","status")]
alluvial_unknown[which(alluvial_unknown$Phylum==""),"Phylum"] <- "Unclassified"
bacteria_alluvial <- alluvial_unknown %>% group_by(Kingdom,Phylum) %>%  count() %>%filter(Kingdom=="Bacteria")%>%arrange(desc(n))
archaea_alluvial <- alluvial_unknown %>% group_by(Kingdom,Phylum) %>%  count() %>%filter(Kingdom=="Archaea")%>%arrange(desc(n))
arc_top <- archaea_alluvial[,"Phylum"]

top10 <- (alluvial_freq %>% group_by(Phylum) %>% summarise(count=sum(n)) %>% arrange(desc(count)) %>% pull(Phylum))[1:10]
source_target1_tmp <- alluvial_unknown %>% mutate(Phylum=case_when((!Phylum %in% top10)~"Others",.default = Phylum)) %>%  group_by(Kingdom,Phylum) %>%  count() %>%arrange(desc(n)) %>% as.data.frame()
colnames(source_target1_tmp) <- c("source","target","value")

source_target2_tmp <- alluvial_unknown %>% mutate(Phylum=case_when((!Phylum %in% top10)~"Others",.default = Phylum)) %>%  group_by(Phylum,status) %>%  count() %>%arrange(desc(n)) %>%as.data.frame()
# source_target2_tmp <- source_target2_tmp %>% 
#   bind_rows(data.frame(Phylum=c("Others","Others"),status=c("Genus+","Species"),n=alluvial_unknown %>%  filter(Kingdom=="Bacteria" & Phylum %in% (bacteria_alluvial[-c(1:15),]%>% pull(Phylum))) %>% group_by(status) %>% count() %>% pull(n)))
# source_target2_tmp <- source_target2_tmp %>% 
#   bind_rows(alluvial_unknown %>% group_by(Phylum,status) %>% filter(Kingdom=="Archaea") %>% count() %>% as.data.frame())

# source_target2_tmp <- alluvial_unknown %>% group_by(Phylum,status) %>%  count() %>%as.data.frame()
colnames(source_target2_tmp) <- c("source","target","value")
source_target_tmp <- source_target1_tmp %>% bind_rows(source_target2_tmp) %>% as.data.frame()
colnames(source_target_tmp) <- c("source","target","value")

nodes_tmp <- data.frame(
  name=c(as.character(source_target_tmp$source), 
         as.character(source_target_tmp$target)) %>% unique()
)
source_target_tmp$IDsource <- match(source_target_tmp$source, nodes_tmp$name)-1 
source_target_tmp$IDtarget <- match(source_target_tmp$target, nodes_tmp$name)-1
library(networkD3)
my_color='d3.scaleOrdinal() .domain(["Archaea","Bacteria","Chloroflexota", "Patescibacteria", "Planctomycetota", "Pseudomonadota", "Nanoarchaeota", "Bacteroidota", "Marinisomatota", "Hydrogenedentota", "Zixibacteria", "Omnitrophota", "Others","Species","Genus+"]).range(["#8b0000","#009900","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#cccccc","#3B5998","#F1C40F"])'
nw_tmp <- sankeyNetwork(Links = source_target_tmp, Nodes = nodes_tmp,
                        Source = "IDsource", Target = "IDtarget",
                        Value = "value", NodeID = "name", 
                        sinksRight=FALSE,
                        width = 500,
                        height = 500,
                        fontSize=12,nodeWidth=10,nodePadding = 6,
                        colourScale = my_color,
                        fontFamily="Arial")
saveNetwork(nw_tmp,"./Prokayotes_Classification_Sankey_tmp_1106_v2.html")

find_lowest_level <- function(x){
  
  ranks <- colnames(x)[1:7]
  #print(ranks)
  x$rank <- NA
  for(l in 2:length(ranks)){
    print(ranks[l])
    x <- mutate(x,rank=case_when((is.na(x[,ranks[l]]) & !is.na(x[ranks[l-1]]))~ ranks[l-1],.default = rank))
  }
  x <- mutate(x,rank=case_when(!is.na(x[,ranks[7]])~ranks[7],.default = rank))
              
  return(x)
}

library(reshape2)
library(ggplot2)

annotation_rank <- summary %>% column_to_rownames(var='user_genome') %>% find_lowest_level(.)

annotation_rank$rank <- factor(annotation_rank$rank,levels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
rc <- annotation_rank%>% group_by(rank) %>% summarise(Freq=n())
rc$Known <-0
rc$Unclassified <- 0
for (taxa in 2:nrow(rc)) {
  rc$Known[taxa] <- nrow(annotation_rank)-sum(rc$Freq[1:taxa-1])
  rc$Unclassified[taxa] <- sum(rc$Freq[1:taxa-1])
}
rc <- rc[-1,-2]
rc$Known <- rc$Known/nrow(annotation_rank)*100
rc$Unclassified <- rc$Unclassified/nrow(annotation_rank)*100
melt_rc <- melt(rc,id.vars = "rank")
melt_rc$rank <- factor(melt_rc$rank,levels = rev(c("Phylum","Class","Order","Family","Genus","Species")))
melt_rc$variable <- factor(melt_rc$variable,levels = c("Unclassified","Known"))
pdf("./MAG_Annotation_Status_20240524.pdf",width = 9,height = 5.5)
ggplot(melt_rc,aes(x=rank,y=value,fill=variable))+geom_bar(stat="identity")+xlab("Taxonomy Ranks")+ylab("Percentage (%)")+
  scale_fill_manual(values = c('#99CCEE',"#99BBEE"))+scale_x_discrete(expand = c(.1,0))+scale_y_continuous(limits =c(-1, 103) ,expand = c(0,0))+
  theme(legend.title = element_blank(),legend.position = "right",legend.text = element_text(size=10),
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size=10),axis.text.y=element_text(size=10,face = "bold"),
        legend.key.size = unit(0.15,"inches"),panel.grid = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black",fill=NA),
        axis.ticks.y   =element_blank() )+
  ggtitle("Total SGBs: 7564")+
  geom_text(aes(label = paste0(round(value,2),"%"),y=value),position=position_stack(vjust = 0.7),size=3,color="#252525")+
  coord_flip()
dev.off()
  