
##packages
library(reshape2)
library(ggClusterNet)
library(tidyverse)
library(patchwork)
library(netET)
library(phyloseq)
#section 3.1
#boxplot for alpha diversity
ggplot(data = PLSO_trait_div1,aes(x=Class,y=values))+
  geom_boxplot(aes(fill=Class),width=0.5,color='black')+
  geom_jitter(position=position_dodge(0.2))+theme_bw()+
  mytheme1+facet_wrap(~index+group,scales = 'free',ncol = 3)+
  scale_fill_manual(values = c("#A6A6A6","#D08441","#7C2F00","#470D00"))+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
#plot for effect size
ggplot(cons_shannon_effects,aes(x=estimate,y=Subgroup))+
  geom_point(aes(color=Subgroup),size=3)+
  geom_errorbarh(aes(xmax = estimate +se, 
                     xmin = estimate-se,
                     color = Subgroup),linewidth=1,height=0.15)+
  geom_text(aes(x=estimate,y=Subgroup,label=sig,hjust=-0.2,vjust=-0.6),size=5)+
  scale_color_manual(values=c( mild= '#D08441',
                               moderate='#7C2F00',
                               severe='#470D00'))+
  geom_vline(xintercept = 0, color = 'gray', linewidth = 1,linetype='dashed')+
  labs(y='',x='Effect Size')+theme_bw()+mytheme1

#PCoA plots
ggplot(PLSO_trait_beta1)+
  geom_point(aes(x=PCoA1,y=PCoA2,color=Class),size=1.5)+
  stat_ellipse(aes(x=PCoA1,y=PCoA2,
                   color=Class,fill=Class,group=Class),
               geom = 'polygon',alpha=0.3,level = 0.9, show.legend = T)+
  scale_color_manual(values=c( "#A6A6A6","#D08441","#7C2F00","#470D00"))+
  scale_fill_manual(values=c( "#A6A6A6","#D08441","#7C2F00","#470D00"))+
  labs(x = 'PCoA1',y = 'PCoA2') +theme_bw()+mytheme1+
  facet_wrap(~group,scales = 'free',ncol = 3)

#section 3.2 co-occurrence network analysis
#Filtering otu by occurrence and minimum abundance
otu_f <- otu|>select_otus(occur_thres = 0.5,abund_thres = 0.001)
#Constructing phyloseq objects
otu_ps <- phyloseq(sample_data(sample_group), #group info
                   otu_table(as.matrix(otu_f), 
                             taxa_are_rows=FALSE), 
                   tax_table(as.matrix(taxonomy)))#taxa table


#Setting the result output path
path='./network_results1/'
#Creating the above path
dir.create(path)
#The network is computed and 
#the resulting edge and node files are used to draw the network graph in Gephi.
network_results <- network.2(ps=otu_ps,
                             maxnode = 4,N=0,
                             big = T,select_layout = TRUE,
                             layout_net = "model_Gephi.2",
                             r.threshold=0.6,p.threshold=0.05,label = FALSE,
                             path = path,zipi = FALSE,
                             group = 'group',ncpus = 4)



##section 3.3
#firstly, mantel test analysis
trait_wood_mantel1 <- 
  mantel_test(t(trait_noc_manteldata),
              wood_properties1[,-1],
              spec_select= list(Consumer=1:2657,
                                Parasite=2658:2867,
                                Phototroph=2868:3297)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
#mantel test-correlation matrix mixplot

trait_wood_mantelplot1 <- 
  qcorrplot(correlate(wood_properties1[,-1],
                      pro_noc_netproperties1[,-1],
                      method = 'spearman'))+geom_tile()+
  geom_mark(size=3,only_mark = T)+
  geom_couple(aes(colour = pd, size = rd),
              data = trait_wood_mantel1, 
              curvature = 0)+scale_size_manual(values = c(0.5, 1, 2))+
  scale_fill_gradientn(colors = c('blue','white','red'),limits=c(-1,1))+
  scale_colour_manual(values = color_pal(3))+
  expand_limits(x=20)+
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"),
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

#section 3.4 Construction of plspm
DC=c(0,0,0,0,0,0,0,0)
WP=c(1,0,0,0,0,0,0,0)
WN=c(1,0,0,0,0,0,0,0)
SP=c(1,1,0,0,0,0,0,0)
SN=c(1,0,1,0,0,0,0,0)
PA=c(1,1,1,1,1,0,0,0)
PC=c(1,1,1,1,1,0,0,0)
PN=c(1,1,1,1,1,0,0,0)
xyc_path7 <- rbind(DC,WP,WN,SP,SN,PA,PC,PN)
xyc_modes7 <- rep('A',8)
cyz2310_blocks7 <- list(1,2,3,4,5,6:7,8:15,21,22:24,25,26)
xyc_plspm_result7 <- plspm(xyc_plspm_testdata7,
                           xyc_path7,xyc_blocks7,
                           modes = xyc_modes7)