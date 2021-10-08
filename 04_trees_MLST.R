library(ggplot2)
library(ggtree) #devtools::install_github(c("hadley/ggplot2", "GuangchuangYu/ggtree"))
library(treeio)
library(ape)
library(scales) #to get ggtree colours
library(dplyr)
library(ggpubr) #ggarrange()
library(ggnewscale)
library(gdsfmt)
library(SNPRelate)
library('dendextend')
library('circlize')
library(ggmap)
library(ggspatial)
library(dplyr)
require(ggrepel)

# inspiration for the code taken from:
# https://yulab-smu.top/treedata-book/
# https://guangchuangyu.github.io/ggtree-book/short-introduction-to-r.html
# https://bioconnector.github.io/workshops/r-ggtree.html

# prevent scientific notation in axis, eg 1E6
options(scipen=10000)

# load MLST
genotype <- read.table("allelic_profiles_ST.txt", stringsAsFactor=F,header = T)
rownames(genotype) <- genotype[,1]
genotype[,1] <- NULL
genotype[] <- lapply( genotype, factor)

# colour categories of sample names
sample_colours <- as.data.frame(read.csv("sample_colors.txt",header = T,sep="\t",stringsAsFactors = F))

# colors to use for MLST
manualcolors<-c('black','dodgerblue2', 'red2', 'gold1', 'cornflowerblue', 
                'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                '#ffff00','firebrick4',  'yellowgreen', 'blue1', 'tan3',
                "maroon4",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
                '#275000', 'moccasin', 'mediumvioletred','cadetblue1',
                "#7e47ff" ,"palegreen4" ,   "tomato3" , "#7CE3D8")


#####PhyML tree

# load tree with bootstrap support
#chr_tree1 <- read.tree("tree_snpsite_267.nwk")
#chr_tree1 <- read.tree("complete_deletion/island_mainland_all_alignment_complete_deletion.phy_phyml_tree_GTR+I+G.txt")
#chr_tree1 <- read.newick("../phyml_267_100bstp_k80/island_mainland_all_trimmed_267.phy_phyml_tree.txt", node.label='support')
chr_tree1 <- read.newick("complete_deletion/island_mainland_all_alignment_complete_deletion.phy_phyml_tree_GTR+I+G.txt", node.label='support')

# create tree:
# change aesthetics of tree (layout, font sizes etc)
# plot root
# plot samples names at the same position
# plot bootstrap support (%)
len_tree <- ggtree(chr_tree1, color="black", size=1.5, linetype=1,  layout = "rectangular", right = F, ladderize = T) +
  theme_tree2(legend.position="none",axis.text.x=element_text(size=16),text=element_text(size=14)) +
  geom_rootedge(rootedge = 0.00005) +
  geom_text(aes(x=branch,label=round(support/10)), size=10*0.36, vjust=-0.5,nudge_x = -0.000022) +
  scale_x_continuous(n.breaks = 6) ;len_tree

# add custom colors to name labels
# highlight our samples
# label clades with samples
len_tree2 <- len_tree %<+% sample_colours + 
  geom_tiplab(size=14*0.36, align=T, linesize=.5,aes(color=factor(my_clr))) +
  scale_color_manual(values = c("black","#0072B2","#D55E00")) +
  geom_hilight(node=70, fill="steelblue", alpha=.6) + 
  geom_hilight(node=75, fill="steelblue", alpha=.6)+ 
  geom_hilight(node=23, fill="steelblue", alpha=.6) +
  geom_text(aes(x=0.00016,y=32,label="group 1"),size=14*0.36) +
  geom_text(aes(x=0.00042,y=30,label="group 2"),size=14*0.36) +
  geom_text(aes(x=0.00029,y=22,label="group 3"),size=14*0.36); len_tree2



# add MLST to the tree, the typing will be added to the correct sample dynamically
pp1 <- (len_tree2) %>%
  gheatmap(genotype, colnames=F, font.size=14*0.36,offset=.0004, width=0.6) + 
  scale_x_ggtree() + 
  scale_fill_manual(name="Housekeeping gene types",values = manualcolors,na.value="grey98") +
  guides(color = "none") ; pp1

# save plot
png("island_mainland_MLST_complete_deletion2.0_phyml.png", res = 300,width =5000,height=3000)
pp1
dev.off()



######## Stricter data ######

# load SNP data
bed.fn <- "../stringent_bed/all_chrom_SNP_stringent_mind2.bed" ; fam.fn <- "../stringent_bed/all_chrom_SNP_stringent_mind2.fam" ; bim.fn <- "../stringent_bed/all_chrom_SNP_stringent_mind2.bim"

# create map object
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "HapMap.gds",cvt.chr="char",cvt.snpid = "auto")
genofile <- snpgdsOpen("HapMap.gds")
snpgdsClose(genofile)

# extract genetic distance matrix
RV <- snpgdsIBS(genofile,snp.id = NULL,verbose = T,autosome.only = F)
m <- 1 - RV$ibs
colnames(m) <- rownames(m) <- RV$sample.id
GeneticDistance <- as.dist(m)

# create hclust object
HC1 <- hclust(GeneticDistance,method="average")

# plot tree
# increase size of plot
# plot sample names
p1 <- ggtree(HC1) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6),legend.position="none") +
  geom_rootedge(rootedge = 0.05) +
  geom_tiplab(size=4, align=T, linesize=.5); p1



pdf("stringent_tree2.pdf", width = 16)
p1
dev.off()



#### BEAST TREES

# 1/x popsize prior

# load beast tree file
beast_tree_cr_norm_rh_none_pops_1x <- read.beast("island_mainland_all_alignment_complete_deletion_cr_norm_rh_none_pops_1x.tree")

# create beast tree with ggtree:
# set most recent sample to 2019
# use custom colours
# plot root
# increase size of plot
# plot sample names
# highlight our samples
# label the clades with our samples
# mark nodes with > 90% posterior support

beast_tree1 <- ggtree(beast_tree_cr_norm_rh_none_pops_1x,mrsd = "2019-01-01", layout = "rectangular", right = F, ladderize = T) ; beast_tree1
beast_tree11 <- beast_tree1 %<+% sample_colours + 
  geom_tiplab(size=14*0.36, align=T, linesize=.5,aes(color=factor(my_clr))) +
  scale_color_manual(values = c("black","#0072B2","#D55E00")) +
  guides(color = "none") +
  geom_rootedge(rootedge = 20) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 140, 6, 6),legend.position="none",text=element_text(size=20)) +
  geom_hilight(node=76, fill="steelblue", alpha=.6) +
  geom_hilight(node=68, fill="steelblue", alpha=.6) +
  geom_hilight(node=32, fill="steelblue", alpha=.6) +
  geom_text(aes(x=1450,y=37,label="group 1"),size=14*0.36) +
  geom_text(aes(x=1650,y=29,label="group 2"),size=14*0.36) +
  geom_text(aes(x=1600,y=21,label="group 3"),size=14*0.36) +
  geom_text(aes(x=-250,y=40,label="A)"),size=14) +
  scale_x_continuous(n.breaks = 10)+ 
  geom_nodepoint(aes(fill=round(as.numeric(posterior)), subset=as.numeric(posterior) > 0.9 ),size=3.5)  ;beast_tree11


# subtract 2019 from heights and multiply by -1 to bring it to CE 
# round up the values
# reverse order of elements, so most recent is second
range_list<-beast_tree11$data$height_0.95_HPD
range_present<-mapply("-",range_list,2019,SIMPLIFY=FALSE)
range_present<-mapply("*",range_present,-1,SIMPLIFY=FALSE)
range_round_ce<-mapply(round,range_present,SIMPLIFY=FALSE)
range_final<-mapply(rev,range_round_ce,SIMPLIFY=FALSE)

# replace in the tree
beast_tree11$data$height_0.95_HPD <- range_final

# set all height values and height ranges that aren't in the nodes we're interested in to NA
beast_tree11$data$height_median[-c(65,66,75,76,67)]<-NA
beast_tree11$data$height_0.95_HPD[-c(65,66,75,76,67)]<-list(c(NULL,NULL))

# and plot height nodes and height ranges
beast_tree111 <- beast_tree11 + geom_nodelab(aes(label=2019-round(height_median, 0)), hjust = -0.05,size=5) +
  geom_nodelab(aes( label= height_0.95_HPD),geom = "text",nudge_y=-1,nudge_x = 40,size=4); beast_tree111

# Lognorm popsize prior


beast_tree_cr_norm_rh_none_pops_lognorm <- read.beast("island_mainland_all_alignment_complete_deletion_cr_norm_rh_none_pops_lognorm.tree")

beast_tree2 <- ggtree(beast_tree_cr_norm_rh_none_pops_lognorm,mrsd = "2019-01-01",  layout = "rectangular", right = F, ladderize = T) ; beast_tree2
beast_tree22 <- beast_tree2 %<+% sample_colours + 
  geom_tiplab(size=14*0.36, align=T, linesize=.5,aes(color=factor(my_clr))) +
  scale_color_manual(values = c("black","#0072B2","#D55E00")) +  guides(color = "none") +
  geom_rootedge(rootedge = 20) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 120, 6, 6),legend.position="none",text=element_text(size=20)) +
  geom_hilight(node=79, fill="steelblue", alpha=.6) +
  geom_hilight(node=67, fill="steelblue", alpha=.6) +
  geom_hilight(node=32, fill="steelblue", alpha=.6) +
  geom_text(aes(x=1750,y=37,label="group 1"),size=14*0.36) +
  geom_text(aes(x=1850,y=28,label="group 2"),size=14*0.36) +
  geom_text(aes(x=1850,y=21,label="group 3"),size=14*0.36) +
  geom_text(aes(x=1000,y=40,label="B)"),size=14) +
  scale_x_continuous(n.breaks = 10) + 
  geom_nodepoint(aes(fill=round(as.numeric(posterior)), subset=as.numeric(posterior) > 0.9 ),size=3.5); beast_tree22


range_list2<-beast_tree22$data$height_0.95_HPD
range_present2<-mapply("-",range_list2,2019,SIMPLIFY=FALSE)
range_present2<-mapply("*",range_present2,-1,SIMPLIFY=FALSE)
range_round_ce2<-mapply(round,range_present2,SIMPLIFY=FALSE)
range_final2<-mapply(rev,range_round_ce2,SIMPLIFY=FALSE)

beast_tree22$data$height_0.95_HPD <- range_final2

beast_tree22$data$height_median[-c(66,65,77,78,67)]<-NA
beast_tree22$data$height_0.95_HPD[-c(66,65,77,78,67)]<-list(c(NULL,NULL))

beast_tree222 <- beast_tree22 + geom_nodelab(aes(label=2019-round(height_median, 0)), hjust = -0.05,size=5) +
  geom_nodelab(aes( label= height_0.95_HPD),geom = "text",nudge_y=-1,nudge_x = 15,size=4); beast_tree222

# combine beast trees and save
beast_trees <- ggarrange(beast_tree111,beast_tree222,ncol = 1)

#png("beast_tree_cr_norm_rh_none_pops_lognorm_+_1x.png", res = 300,width =5000,height=6500)
pdf("beast_tree_cr_norm_rh_none_pops_lognorm_+_1x_group.pdf",height = 20,width=15)
beast_trees
dev.off()




### SAMPLE LOCATIONS

# set up where the western isles are on the map
bounds <- c(left =  -8, 
            bottom =57,
            right = -5, 
            top = 58.5)

# load site sample information metadata
samples<-read.csv("site_locations.txt",header = T,sep = "\t")

# create terrain map
map <- get_stamenmap(bounds,zoom=10, maptype = "terrain", invert = T)

# plot map
# plot site locations
# change color of points; 1 for north one for south; saved in file
# add labels with location names
# remove axis and legend
# add title 
# add north notation
location_plot <- ggmap(map) + 
  geom_point(data = samples, 
    mapping = aes(x = lon, y = lat, fill = colour), size=12,shape=21) +
  
  geom_text_repel(data = samples,aes(x = lon, y = lat, label = site), 
    size = 7.5,box.padding = 0.43) +
  
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, -1, -1), 'lines'),
        legend.position="none",
        text=element_text(size=20)) +
  
  xlab('') +
  ylab('') +
  
  ggtitle("Sampling locations in the Western Isles") +
  
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.3, "in"), pad_y = unit(8.5, "in"), 
                         style = north_arrow_fancy_orienteering); location_plot


png("location_map.png", res = 300,width =5000,height=3000)
location_plot
dev.off()


