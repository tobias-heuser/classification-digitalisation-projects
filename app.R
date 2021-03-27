# install.packages("ggdendro")
# install.packages("reshape2")
# install.packages("grid")
# install.packages("dendextend")

library(dplyr)

projects <- data.frame(
  #id.s = c(1:58),
  pro_involved = c("+C", "1D", "+D", "+C", "+C", "+D", "+D", "+D", "1D", "1D", "+C", "+D", "+C", "+D", "+C", "+C", "+C", "+D", "+C", "+C", "+C", "1D", "+D", "+C", "1D", "+C", "+D", "+C", "+C", "+C", "+C", "+D", "+C", "+D", "+D", "+D", "+D", "+C", "+D", "+C", "+C", "+C", "+C", "+C", "+C", "+D", "+C", "+D", "+C", "+C", "+C", "+C", "+C", "+C", "+C", "+C", "+C", "+C"),
  pro_focus = c("P&S", "OE&A", "I&DM", "I&DM", "C", "OE&A", "OE&A", "OE&A", "OE&A", "OE&A", "P&S", "P&S", "C", "OE&A", "P&S", "OE&A", "C", "I&DM", "C", "I&DM", "OE&A", "I&DM", "C", "P&S", "I&DM", "C", "OE&A", "OE&A", "I&DM", "I&DM", "OE&A", "I&DM", "I&DM", "OE&A", "OE&A", "OE&A", "OE&A", "C", "I&DM", "I&DM", "I&DM", "P&S", "I&DM", "P&S", "I&DM", "OE&A", "C", "P&S", "P&S", "OE&A", "OE&A", "C", "OE&A", "I&DM", "OE&A", "OE&A", "I&DM", "P&S"),
  pro_complexity = c("H/H", "L/L", "L/H", "H/H", "H/H", "H/L", "H/L", "H/H", "H/H", "H/H", "L/H", "H/L", "L/H", "H/H", "L/H", "H/H", "L/H", "H/H", "H/H", "H/H", "L/H", "H/H", "H/L", "H/H", "H/L", "H/H", "H/H", "H/L", "H/L", "H/H", "L/L", "L/H", "H/L", "L/L", "H/H", "L/L", "H/H", "H/H", "L/H", "L/H", "H/H", "H/H", "H/H", "L/H", "H/H", "H/H", "H/L", "H/L", "H/H", "H/L", "L/L", "L/L", "H/L", "L/L", "H/L", "H/L", "L/H", "H/H"), 
  pro_impact = c("+R", "-C", "+R", "-C", "+R", "-C", "+R", "-C", "O", "O", "+R", "O", "+R", "O", "O", "-C", "+R", "O", "+R", "O", "-C", "O", "+R", "+R", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "-C", "-C", "-C", "+R", "O", "-C", "O", "+R", "O", "O", "O", "-C", "+R", "+R", "+R", "+R", "O", "-C", "-C", "-C", "-C", "-C", "O", "+R"), 
  pro_mode = c("explore", "explore", "explore", "exploit", "exploit", "exploit", "explore", "exploit", "exploit", "exploit", "explore", "exploit", "exploit", "exploit", "exploit", "exploit", "exploit", "exploit", "explore", "explore", "exploit", "exploit", "explore", "explore", "exploit", "explore", "exploit", "exploit", "exploit", "explore", "explore", "exploit", "exploit", "exploit", "exploit", "exploit", "exploit", "explore", "exploit", "exploit", "explore", "exploit", "exploit", "exploit", "explore", "exploit", "explore", "explore", "exploit", "exploit", "exploit", "exploit", "exploit", "exploit", "exploit", "exploit", "exploit", "explore"),
  pro_pathway = c("CX", "CX", "I&S", "I&S", "CX", "iterate", "new", "I&S", "I&S", "I&S", "CX", "CX", "CX", "I&S", "new", "new", "new", "I&S", "new", "I&S", "I&S", "I&S", "CX", "I&S", "I&S", "new", "I&S", "I&S", "new", "new", "I&S", "I&S", "I&S", "I&S", "I&S", "I&S", "I&S", "CX", "CX", "I&S", "CX", "CX", "new", "CX", "CX", "I&S", "CX", "CX", "CX", "I&S", "CX", "I&S", "new", "I&S", "I&S", "I&S", "I&S", "iterate"),
  stringsAsFactors=TRUE
)

#----- Dissimilarity Matrix -----#
library(cluster) 
# to perform different types of hierarchical clustering
# package functions used: daisy(), diana(), clusplot()
gower.dist <- daisy(projects[ ,1:6], metric = c("gower"))
# class(gower.dist) 
## dissimilarity , dist

#------------ DIVISIVE CLUSTERING ------------#
divisive.clust <- diana(as.matrix(gower.dist), 
                        diss = TRUE, keep.diss = TRUE)
plot(divisive.clust,
     main = "Divisive")

#------------ AGGLOMERATIVE CLUSTERING ------------#
aggl.clust.c <- hclust(gower.dist, method = "complete")
plot(aggl.clust.c,
     main = "Agglomerative, complete linkages")

# Cluster stats comes in a list form, it is more convenient to look at it as a table
# This code below will produce a dataframe with observations in columns and variables in row
# Not quite tidy data, but it's nicer to look at
library(fpc)

cstats.table <- function(dist, tree, k) {
  clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between",
                    "wb.ratio","dunn2","avg.silwidth")
  clust.size <- c("cluster.size")
  stats.names <- c()
  row.clust <- c()
  
  output.stats <- matrix(ncol = k, nrow = length(clust.assess))
  cluster.sizes <- matrix(ncol = k, nrow = k)
  
  for(i in c(1:k)){
    row.clust[i] <- paste("Cluster-", i, " size")
  }
  
  for(i in c(2:k)){
    stats.names[i] <- paste("Test", i-1)
    
    for(j in seq_along(clust.assess)){
      output.stats[j, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.assess])[j]
      
    }
    
    for(d in 1:k) {
      cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
      dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
      cluster.sizes[d, i]
      
    }
  }
  
  output.stats.df <- data.frame(output.stats)
  
  cluster.sizes <- data.frame(cluster.sizes)
  cluster.sizes[is.na(cluster.sizes)] <- 0
  
  rows.all <- c(clust.assess, row.clust)
  # rownames(output.stats.df) <- clust.assess
  output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
  colnames(output) <- stats.names[2:k]
  rownames(output) <- rows.all
  
  is.num <- sapply(output, is.numeric)
  output[is.num] <- lapply(output[is.num], round, 2)
  
  output
}

# I am capping the maximum amout of clusters by 7
# but for sure, we can do more
# I want to choose a reasonable number, based on which I will be able to see basic differences between customer groups

stats.df.divisive <- cstats.table(gower.dist, divisive.clust, 10)
stats.df.divisive

stats.df.aggl <- cstats.table(gower.dist, aggl.clust.c, 10)
stats.df.aggl

# --------- Choosing the number of clusters ---------#

# Using "Elbow" and "Silhouette" methods to identify the best number of clusters

library(ggplot2)

# Elbow
# Divisive clustering
ggplot(data = data.frame(t(stats.df.divisive)), aes(x=cluster.number, y=within.cluster.ss)) + geom_point()+
  geom_line()+
  ggtitle("") +
  labs(x = "Num.of clusters", y = "Within sum of squares") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw(base_size=20)

# Silhouette
ggplot(data = data.frame(t(stats.df.divisive)), aes(x=cluster.number, y=avg.silwidth)) + geom_point()+
  geom_line()+
  ggtitle("Divisive clustering") +
  labs(x = "Num.of clusters", y = "Average silhouette width") +
  theme(plot.title = element_text(hjust = 0.5))

# Agglomorative clustering
# Elbow
ggplot(data = data.frame(t(stats.df.aggl)), aes(x=cluster.number, y=within.cluster.ss)) + geom_point()+
  geom_line()+
  ggtitle("") +
  labs(x = "Num.of clusters", y = "Within sum of squares") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw(base_size=20)

# Silhouette
ggplot(data = data.frame(t(stats.df.aggl)), aes(x=cluster.number, y=avg.silwidth)) + geom_point()+
  geom_line()+
  ggtitle("Agglomorative clustering") +
  labs(x = "Num.of clusters", y = "Average silhouette width") +
  theme(plot.title = element_text(hjust = 0.5))

# Finally, assigning the cluster number to the observation

clust.num <- cutree(divisive.clust, k = 3) 
id.s = c(1:58)
projects.cl <- cbind(id.s, projects, clust.num)

clust.aggl.num <- cutree(aggl.clust.c, k = 3) 
id.s = c(1:58)
projects.aggl.cl <- cbind(id.s, projects, clust.aggl.num)
#projects.cl <- cbind(projects, clust.num)

library("ggplot2")
library("reshape2")
library("purrr")
library("dplyr")
library("dendextend")

dendro <- as.dendrogram(aggl.clust.c)

dendro.col <- dendro %>%
  set("branches_k_color", k = 3, 
      value = c("gold3", "darkcyan", "cyan3")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors", 
      value = c("darkslategray")) %>% 
  set("labels_cex", 0.5) 


ggd1 <- as.ggdend(dendro.col)
ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram (aggl), k = 3")


# Create a radial plot
ggplot(ggd1, labels = T) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")

# cust.order <- order.dendrogram(dendro)
# projects.cl.ord <- projects.cl[cust.order, ]

# 1 variable per row
# factors have to be converted to characters in order not to be dropped

cust.long <- melt(data.frame(lapply(projects.cl, as.character), stringsAsFactors=FALSE), 
                  id.vars = c("id.s", "clust.num"), factorsAsStrings=T)

cust.aggl.long <- melt(data.frame(lapply(projects.aggl.cl, as.character), stringsAsFactors=FALSE), 
                  id.vars = c("id.s", "clust.aggl.num"), factorsAsStrings=T)

cust.long.q <- cust.long %>%
  group_by(clust.num, variable, value) %>%
  mutate(count = n_distinct(id.s)) %>%
  distinct(clust.num, variable, value, count)

cust.aggl.long.q <- cust.aggl.long %>%
  group_by(clust.aggl.num, variable, value) %>%
  mutate(count = n_distinct(id.s)) %>%
  distinct(clust.aggl.num, variable, value, count)

cust.long.p <- cust.long.q %>%
  group_by(clust.num, variable) %>%
  mutate(perc = count / sum(count)) %>%
  arrange(clust.num)

cust.aggl.long.p <- cust.aggl.long.q %>%
  group_by(clust.aggl.num, variable) %>%
  mutate(perc = count / sum(count)) %>%
  arrange(clust.aggl.num)

heatmap.p <- ggplot(cust.long.p, aes(x = clust.num, y = factor(value, levels = c("1D","+D","+C",
                                                                                 "C", "P&S", "I&DM", "OE&A",
                                                                                 "L/L","L/H", "H/L", "H/H",
                                                                                 "+R","-C","O",
                                                                                 "exploit","explore",
                                                                                 "I&S","CX","iterate","new"), 
                                                               ordered = T))) +
  geom_tile(aes(fill = perc), alpha = 0.85)+
  labs(title = "Distribution of characteristics across clusters", x = "Cluster number", y = NULL) +
  geom_hline(yintercept = 3.5) + 
  geom_hline(yintercept = 7.5) + 
  geom_hline(yintercept = 11.5) + 
  geom_hline(yintercept = 14.5) + 
  geom_hline(yintercept = 16.5) + 
  scale_fill_gradient2(low = "darkslategray1", mid = "yellow", high = "turquoise4")

heatmap.p

heatmap.aggl.p <- ggplot(cust.aggl.long.p, aes(x = clust.aggl.num, y = factor(value, levels = c("1D","+D","+C",
                                                                                 "C", "P&S", "I&DM", "OE&A",
                                                                                 "L/L","L/H", "H/L", "H/H",
                                                                                 "+R","-C","O",
                                                                                 "exploit","explore",
                                                                                 "I&S","CX","iterate","new"), 
                                                               ordered = T))) +
  geom_tile(aes(fill = perc), alpha = 0.85)+
  labs(title = "Distribution of characteristics across (aggl) clusters", x = "Cluster number", y = NULL) +
  geom_hline(yintercept = 3.5) + 
  geom_hline(yintercept = 7.5) + 
  geom_hline(yintercept = 11.5) + 
  geom_hline(yintercept = 14.5) + 
  geom_hline(yintercept = 16.5) + 
  scale_fill_gradient2(low = "darkslategray1", mid = "yellow", high = "turquoise4")

heatmap.aggl.p
