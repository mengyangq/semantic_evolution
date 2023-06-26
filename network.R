library(coop)
library(NetworkToolbox)
library(igraph)
library(dplyr)
library(RColorBrewer)
library(Polychrome)

# read the binary matrix for 'stone'
stone <- read.csv("data/stone.csv")

# read the detailed info for each concept/node in 'stone' (some of the info were appended after the network has been created, e.g., degree)
stone_map <- read.csv("data/stone_map.csv")

# calculate cosine similarity for each pair of concpets in 'stone'
stone_cosine <- cosine(stone)

# TMFG filtering
stone_tmfg <- TMFG(stone_cosine, normal = FALSE)$A

# create the network
stone_graph <- graph_from_adjacency_matrix(stone_tmfg != 0, mode = "undirected", diag = FALSE)


# obtain the degree of each node and create the degree distribution plot
degree_stone <- degree(stone_graph)
degree_stone <- table(degree_stone)

relafreq_stone <- degree_stone / sum(degree_stone)

barplot(relafreq_stone,
  xlab = "Degree", ylab = "Relative Frequencies",
  main = "Stone",
  col = "gray60",
  ylim = range(pretty(c(0, relafreq_stone)))
)

# obtain the eigenvector centrality
write.csv(centr_eigen(stone_graph),'stone_eigen.csv')

# get CC and ASPL
stone_cc <- transitivity(stone_graph, type = "global", vids = NULL, weights = NULL)
stone_mean_distance <- mean_distance(stone_graph, unconnected = TRUE)

# CC and ASPL of 1000 random network simulations
stone_random_cc <- numeric(1000)
stone_random_mean_distance <- numeric(1000)
for (i in 1:1000) {
  rn <- erdos.renyi.game(vcount(stone_graph), ecount(stone_graph), type = "gnm")
  stone_random_cc[i] <- transitivity(rn, type = "global", vids = NULL, weights = NULL)
  stone_random_mean_distance[i] <- mean_distance(rn, unconnected = TRUE)
}

# z-tests of CC and ASPL between stone network and random networks
stone_cc_test <- z.test(stone_random_cc, sigma.x = sd(stone_random_cc), mu = stone_cc)
stone_mean_distance_test <- z.test(stone_random_mean_distance, sigma.x = sd(stone_random_mean_distance), mu = stone_mean_distance)


# change the graph/network back to data frame (vertices and edges) in order to append the additional info to the vertices, then recreate the same network
stone_df <- igraph::as_data_frame(stone_graph, what = "both")

stone_df$vertices <- stone_df$vertices %>%
  left_join(stone_map, c("name" = "Code"))

updated_stone_graph <- graph_from_data_frame(stone_df$edges, directed = F, vertices = stone_df$vertices)

# optimal community structure
stone_cluster <- cluster_optimal(updated_stone_graph)
# write.csv(stone_cluster$membership,'stone_cluster.csv')

# network layout
stone_layout <- layout_(updated_stone_graph, with_graphopt())

# number of unique categories
stone_num_categories <- length(unique(stone_map$Category))

# generate a palette of colors with the 'dark.colors' function
stone_category_colors <- dark.colors(stone_num_categories)

# assign categories to colors
named_stone_category_colors <- setNames(stone_category_colors, sort(unique(stone_map$Category)))
stone_category_indices <- match(stone_map$Category, names(named_stone_category_colors))


# number of unique clusters
stone_num_clusters <- length(unique(stone_map$cluster_membership))

# create a palette with contrasting colors
stone_cluster_colors <- brewer.pal(n = min(stone_num_clusters, 12), name = "Set3")

# assign clusters to colors
named_stone_cluster_colors <- setNames(stone_cluster_colors, sort(unique(stone_map$cluster_membership)))
stone_cluster_indices <- match(stone_map$cluster_membership, names(named_stone_cluster_colors))


# scale the 'degree' column to a range from 1 to 10
degree_scaled <- 1 + 9 * (stone_map$degree - min(stone_map$degree)) / (max(stone_map$degree) - min(stone_map$degree))
degree_ordered <- degree_scaled[order(as.numeric(names(V(updated_stone_graph))))]

# create the network visualization
pdf(file = "stone.pdf", width = 16, height = 14)

plot(updated_stone_graph, vertex.label = V(updated_stone_graph)$Abbreviated, vertex.color = named_stone_cluster_colors[stone_cluster_indices], vertex.label.color = named_stone_category_colors[stone_category_indices], vertex.label.cex = 1.25, vertex.size = degree_ordered, vertex.frame.color = NA, vertex.label.dist = 0.4, edge.width = .1, edge.curved = .1, layout = stone_layout)

# add a title
title("Stone", cex.main = 2)

# add a legend
legend(x = "topleft", legend = names(named_stone_category_colors), fill = named_stone_category_colors, cex = 1.5, title = "Categories", box.lty = 0)
legend(x = "bottomright", legend = names(named_stone_cluster_colors), fill = named_stone_cluster_colors, cex = 1.5, title = "Clusters", box.lty = 0)

dev.off()



# same steps for other networks
