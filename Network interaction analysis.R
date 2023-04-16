library(phyloseq)
library(SpiecEasi)

path <- "/home/simonw/Schreibtisch/Sequencing_Data/Fungiome/Crohn_Kohorte"

setwd("/home/simonw/Schreibtisch/Sequencing_Data/Fungiome/Crohn_Kohorte/")

count_tab= read.csv(file = 'feature-table-crohn_excluded2.csv', sep = '\t', header=TRUE, row.names=1)

head(count_tab)


tax_tab= read.csv(file = 'taxonomy2.csv', sep = '\t', header=TRUE, row.names=1)

head(tax_tab)





tax_tab =as.matrix(tax_tab)
count_tab =as.matrix(count_tab)


physeq1 <- phyloseq(otu_table(count_tab, taxa_are_rows = T),   #taxa_are_rows=F (if your taxa names on the column not the rows)
                   tax_table(tax_tab))


countperphyla = 3
Samplepercentage = 0.05

physeq1 = filter_taxa(physeq1, function(x) sum(x > countperphyla) > (Samplepercentage*length(x)), TRUE)

#### Normalize number of reads in each sample using median sequencing depth.####

total = median(sample_sums(physeq1))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq_mednorm = transform_sample_counts(physeq1, standf)

# Transform to relative abundance. Save as new object.
physeq_re = transform_sample_counts(physeq_mednorm, function(x){x / sum(x)})

physeq_re <- subset_taxa(physeq_re, !is.na(Genus) & !Genus%in% c("", " g__","uncharacterized"))

physeq_re <- subset_taxa(physeq_re, !is.na(Species) & !Species%in% c("", " s__","uncharacterized"))



#se.mb.amgut2 <- spiec.easi(physeq_re, method='mb', lambda.min.ratio=1e-2,
                           #nlambda=20, pulsar.params=list(rep.num=50))
#ig2.mb <- adj2igraph(getRefit(se.mb.amgut2),  vertex.attr=list(name=taxa_names(physeq_re)))
#plot_network(ig2.mb, physeq_re, type='taxa', color="Genus")



#Mycobiome
path <- "/home/simonw/Schreibtisch/Sequencing_Data/Fungiome/Crohn_Kohorte"

setwd("/home/simonw/Schreibtisch/Sequencing_Data/Fungiome/Crohn_Kohorte/")
count_tab= read.csv(file = 'asv_count_fungi_AA83ex_2.csv', sep = '\t', header=TRUE, row.names=1)

head(count_tab)


tax_tab= read.csv(file = 'asv_tax_unite.tsv', sep = '\t', header=TRUE, row.names=1)

head(tax_tab)




tax_tab =as.matrix(tax_tab)
count_tab =as.matrix(count_tab)

physeq2 <- phyloseq(otu_table(count_tab, taxa_are_rows = T),   #taxa_are_rows=F (if your taxa names on the column not the rows)
                    tax_table(tax_tab))



countperphyla = 3
Samplepercentage = 0.05

physeq2 = filter_taxa(physeq2, function(x) sum(x > countperphyla) > (Samplepercentage*length(x)), TRUE)

#### Normalize number of reads in each sample using median sequencing depth.####

total = median(sample_sums(physeq2))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq_mednorm2 = transform_sample_counts(physeq2, standf)

# Transform to relative abundance. Save as new object.
physeq_re2 = transform_sample_counts(physeq_mednorm2, function(x){x / sum(x)})

physeq_re2 <- subset_taxa(physeq_re2, !is.na(Genus) & !Genus%in% c("", " g__","uncharacterized"))

physeq_re2 <- subset_taxa(physeq_re2, !is.na(Species) & !Species%in% c("", " s__","uncharacterized"))



#se.mb.amgut2 <- spiec.easi(physeq_re2, method='mb', lambda.min.ratio=1e-2,
 #                          nlambda=20, pulsar.params=list(rep.num=50))
#ig2.mb <- adj2igraph(getRefit(se.mb.amgut2),  vertex.attr=list(name=taxa_names(physeq_re2)))
#plot_network(ig2.mb, physeq_re2, type='taxa', color="Species")


## Default settings ##

## Parallel multicore ##
pargs2 <- list(rep.num=50, seed=10010, ncores=4)

t4 <- system.time(
  se4 <- spiec.easi(list(physeq_re, physeq_re2), method='glasso', lambda.min.ratio=1e-3, nlambda=30,
                    sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs2)
)

se.hmp2 <- se4 
#se.hmp2 <- spiec.easi(list(physeq_re, physeq_re2), method='glasso', nlambda=40,
                      #lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))

dtype <- c(rep(1,ntaxa(physeq_re)), rep(2,ntaxa(physeq_re2)))
dtype <- c(rep(1,ntaxa(hmp216S)), rep(2,ntaxa(hmp2prot)))
print(dtype)

names <- append(taxa_names(physeq_re2),taxa_names(physeq_re))
names
physeq_all<- merge_phyloseq(physeq_re,physeq_re2)
se1 <- adj2igraph(getRefit(se4),  vertex.attr=list(name=names))
print(se1)


Isolated = which(degree(se1)<=48)
se2 = delete.vertices(se1, Isolated)
plot(se2, layout=LO2)
plot_network(se2, physeq_all, type='taxa', color="Genus",point_size=18,line_alpha=0.5,line_weight=0.05, label="Species",hjust=0.5,alpha=0.1)

library("RCy3")
createNetworkFromIgraph(se1 ,"myIgraph")
write.table(igraph::as_edgelist(se4 , names = T),"Edges.txt")



#Color depending on degree
library(igraph)

grph.pos <- se1

dd.grph.pos <- degree.distribution(grph.pos)
plot(0:(length(dd.grph.pos)-1), dd.grph.pos, type='b',
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

grph.pos_deg<-degree(grph.pos, v=V(grph.pos), mode="all")

fine = 500 # this will adjust the resolving power.

#this gives you the colors you want for every point
library(viridis)
graphCol = viridis(fine)[as.numeric(cut(grph.pos_deg,breaks = fine))]

# now plot

plot.igraph(grph.pos)
plot(grph.pos, vertex.color=graphCol,
     edge.color="black",
     vertex.attr=list(name=names),
     vertex.label=NA,
     vertex.color=dtype+1,
     vertex.size=2.5,
     layout=layout_with_fr(grph.pos))




Sgrph.pos_bw<-betweenness(grph.pos, directed=F)

#this gives you the colors you want for every point
graphCol = viridis(fine)[as.numeric(cut(grph.pos_bw,breaks = fine))]

# now plot
plot(grph.pos, vertex.color=graphCol,
     vertex.label=NA,
     edge.color="black",
     vertex.color=dtype+1,
     vertex.size=3,
     layout=layout_with_fr(grph.pos))

grph.pos_tran<-transitivity(grph.pos, type="local")
grph.pos_tran

#this gives you the colors you want for every point
graphCol = viridis(fine)[as.numeric(cut(grph.pos_tran,breaks = fine))]

# now plot
plot(grph.pos, vertex.color=graphCol,
     vertex.label=NA,
     edge.color="black",
     vertex.color=dtype+1,
     vertex.size=3,
     layout=layout_with_fr(grph.pos))

grph.pos_tran_gl<-transitivity(grph.pos, type="global")
grph.pos_tran_gl


grph.pos.greedy <- cluster_fast_greedy(grph.pos, weights=E(grph.pos)$weight)
modularity(grph.pos.greedy)
sizes(grph.pos.greedy)

colourCount = length(unique(grph.pos.greedy$membership)) # this will adjust the resolving power.

cluster_col = rainbow(colourCount)[as.numeric(cut(grph.pos.greedy$membership,breaks = colourCount))]

plot(grph.pos, vertex.color=cluster_col,
     vertex.label=NA,
     edge.color="black",
     vertex.color=dtype+1,
     vertex.size=2.5,
     layout=layout_with_fr(grph.pos))

grph.pos.louvain <- cluster_louvain(grph.pos, weights=E(grph.pos)$weight)
modularity(grph.pos.louvain)

sizes(grph.pos.louvain)

colourCount = length(unique(grph.pos.louvain$membership)) # this will adjust the resolving power.

cluster_col = rainbow(colourCount)[as.numeric(cut(grph.pos.louvain$membership,breaks = colourCount))]

#plot <-  adj2igraph(getRefit(grph.pos), vertex.color=cluster_col)

plot(grph.pos, vertex.color=cluster_col,
     vertex.label=NA,
     edge.color="black",
     vertex.size=3,
     layout=layout_with_fr(grph.pos))


se.cor  <- cov2cor(as.matrix(getOptCov(se.hmp2)))
weighted.adj.mat <- se.cor*getRefit(se.hmp2)

grph <- adj2igraph(weighted.adj.mat)

grph_whole <- adj2igraph(weighted.adj.mat)
grph_whole<-delete.edges(grph_whole,which(E(grph_whole)$weight<0))
grph.whole.louvain <- cluster_louvain(grph_whole, weights=E(grph_whole)$weight)
modularity(grph.whole.louvain)

sizes(grph.whole.louvain)

colourCount = length(unique(grph.whole.louvain$membership)) # this will adjust the resolving power.

cluster_col = rainbow(colourCount)[as.numeric(cut(grph.whole.louvain$membership,breaks = colourCount))]

# now plot
plot(grph_whole, vertex.color=cluster_col,
     vertex.label=NA,
     edge.color="black",
     vertex.size=4,
     layout=layout_with_fr(grph_whole))

V(grph.pos)$cluster=grph.pos.louvain$membership
vertex_attr(grph.pos, index = V(grph.pos))

ids <- which(sizes(grph.pos.louvain)<=2)
grph.pos.main.communities <- delete_vertices(grph.pos,which(V(grph.pos)$cluster %in% ids))

nodes <- V(grph.pos.main.communities)$name
nodes

write.table(nodes,"Cluster2.csv")

cluster_id <- V(grph.pos.main.communities)$cluster


nodes<-as.data.frame(cbind(nodes, cluster_id))
print(nodes)

library(dplyr)

taxa= read.csv(file = 'taxonomy2_combined.csv', sep = '\t',header=TRUE)
taxa=as.data.frame(taxa)
print(taxa)

colnames(nodes)<-c("Names","Louvain Cluster")
nodes=as.data.frame(nodes)
nodes<-merge(nodes, taxa,  by="Names",all.x=TRUE)

vertex_attr(grph.pos, index = V(grph.pos)) 

nodes

Family_breakdown<-table(nodes$Family,nodes$`Louvain Cluster`) 
Family_breakdown<-as.data.frame(Family_breakdown)

library("ggplot2")

ggplot(Family_breakdown) +
  geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = 'identity', width = 0.5) +
  labs(x = "Louvain cluster",
       y = "Count") +
  guides(fill=guide_legend(title="Louvain cluster")) +
  #scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                               #"#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#661100", "#44AA99")) +
  theme_bw() +
  theme(panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        axis.text.x = element_text(angle=90, hjust=1))
