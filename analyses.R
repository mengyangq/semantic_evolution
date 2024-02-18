library(coop)
library(NetworkToolbox)
library(igraph)
library(dplyr)
library(RColorBrewer)
library(Polychrome)
library(car)
library(BSDA)
library(vcd)

### NETWORK CONSTRUCTION

# read the binary matrix for 'fruit'
fruit <- read.csv("fruit.csv")

# calculate cosine similarity for each pair of concpets in 'fruit'
fruit_cosine <- cosine(fruit)

# TMFG filtering
fruit_tmfg <- TMFG(fruit_cosine, normal = FALSE)$A

# create the network
fruit_graph <- graph_from_adjacency_matrix(fruit_tmfg != 0, mode = "undirected", diag = FALSE)


# obtain the degree of each node and create the degree distribution plot
degree_fruit <- degree(fruit_graph)
degree_fruit <- table(degree_fruit)

relafreq_fruit <- degree_fruit / sum(degree_fruit)

barplot(relafreq_fruit,
        xlab = "Degree", ylab = "Relative Frequencies",
        main = "Fruit",
        col = "gray60",
        ylim = range(pretty(c(0, relafreq_fruit)))
)

write.csv(degree(fruit_graph),'fruit_degree.csv')

# optimal community structure
fruit_cluster <- cluster_optimal(fruit_graph)

# Create a data frame with node IDs and memberships
node_ids <- V(fruit_graph)$name  # assuming nodes are named; if not, use 1:vcount(fruit_graph)
memberships <- fruit_cluster$membership

# Combine into a data frame
cluster_data <- data.frame(Node = node_ids, Membership = memberships)

# Write to CSV
write.csv(cluster_data, 'fruit_cluster.csv', row.names = FALSE)

# get CC and ASPL
fruit_cc <- transitivity(fruit_graph, type = "global", vids = NULL, weights = NULL)
fruit_mean_distance <- mean_distance(fruit_graph, unconnected = TRUE)

# CC and ASPL of 1000 random network simulations
fruit_random_cc <- numeric(1000)
fruit_random_mean_distance <- numeric(1000)
for (i in 1:1000) {
  rn <- erdos.renyi.game(vcount(fruit_graph), ecount(fruit_graph), type = "gnm")
  fruit_random_cc[i] <- transitivity(rn, type = "global", vids = NULL, weights = NULL)
  fruit_random_mean_distance[i] <- mean_distance(rn, unconnected = TRUE)
}

# z-tests of CC and ASPL between fruit network and random networks
fruit_cc_test <- z.test(fruit_random_cc, sigma.x = sd(fruit_random_cc), mu = fruit_cc)
fruit_mean_distance_test <- z.test(fruit_random_mean_distance, sigma.x = sd(fruit_random_mean_distance), mu = fruit_mean_distance)


# read the detailed info for each concept/node in 'fruit' (some of the info were appended after the network has been created, e.g., degree)
fruit_map <- read.csv("fruit_map.csv")

# change the graph/network back to data frame (vertices and edges) in order to append the additional info to the vertices, then recreate the same network
fruit_df <- igraph::as_data_frame(fruit_graph, what = "both")

fruit_df$vertices <- fruit_df$vertices %>%
  left_join(fruit_map, c("name" = "Code"))

updated_fruit_graph <- graph_from_data_frame(fruit_df$edges, directed = F, vertices = fruit_df$vertices)

# network layout
fruit_layout <- layout_(updated_fruit_graph, with_graphopt())

# number of unique categories
fruit_num_categories <- length(unique(fruit_map$Category))

# generate a palette of colors with the 'dark.colors' function
fruit_category_colors <- dark.colors(fruit_num_categories)

# assign categories to colors
named_fruit_category_colors <- setNames(fruit_category_colors, sort(unique(fruit_map$Category)))
fruit_category_indices <- match(fruit_map$Category, names(named_fruit_category_colors))


# number of unique clusters
fruit_num_clusters <- length(unique(fruit_map$Cluster))

# create a palette with contrasting colors
fruit_cluster_colors <- brewer.pal(n = min(fruit_num_clusters, 12), name = "Set3")

# assign clusters to colors
named_fruit_cluster_colors <- setNames(fruit_cluster_colors, sort(unique(fruit_map$Cluster)))
fruit_cluster_indices <- match(fruit_map$Cluster, names(named_fruit_cluster_colors))


# scale the 'degree' column to a range from 1 to 10
degree_scaled <- 1 + 9 * (fruit_map$Degree - min(fruit_map$Degree)) / (max(fruit_map$Degree) - min(fruit_map$Degree))
degree_ordered <- degree_scaled[order(as.numeric(names(V(updated_fruit_graph))))]

# create the network visualization
pdf(file = "fruit.pdf", width = 16, height = 14)

plot(updated_fruit_graph, vertex.label = V(updated_fruit_graph)$Abbreviated, vertex.color = named_fruit_cluster_colors[fruit_cluster_indices], vertex.label.color = named_fruit_category_colors[fruit_category_indices], vertex.label.cex = 1.25, vertex.size = degree_ordered, vertex.frame.color = NA, vertex.label.dist = 0.4, edge.width = .1, edge.curved = .1, layout = fruit_layout)

# add a title
title("Fruit", cex.main = 2)

# add a legend
legend(x = "topleft", legend = names(named_fruit_category_colors), fill = named_fruit_category_colors, cex = 1.5, title = "Categories", box.lty = 0)
legend(x = "bottomright", legend = names(named_fruit_cluster_colors), fill = named_fruit_cluster_colors, cex = 1.5, title = "Clusters", box.lty = 0)

dev.off()

# Similar steps for other networks: stone, languages (fruit), languages (stone) 

### STABILITY

bootstrap_network_metrics <- function(data, n_iterations) {
  cc_results <- numeric(n_iterations)
  aspl_results <- numeric(n_iterations)
  
  for (i in 1:n_iterations) {
    # Bootstrap with replacement
    bootstrapped_sample <- data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
    
    # Calculate cosine similarity and create network
    cosine_matrix <- cosine(bootstrapped_sample)
    cosine_matrix[is.na(cosine_matrix)] <- 0
    tmfg_network <- TMFG(cosine_matrix, normal = FALSE)$A
    graph <- graph_from_adjacency_matrix(tmfg_network != 0, mode = "undirected", diag = FALSE)
    
    # Calculate metrics
    cc_results[i] <- transitivity(graph, type = "global")
    aspl_results[i] <- mean_distance(graph, unconnected = TRUE)

  }
  
  return(list(CC = cc_results, ASPL = aspl_results))
  
}

set.seed(123) # for reproducibility
bootstrap_fruit <- bootstrap_network_metrics(fruit, 1000)

mean_cc_fruit <- mean(bootstrap_fruit$CC)
sd_cc_fruit <- sd(bootstrap_fruit$CC)
mean_aspl_fruit <- mean(bootstrap_fruit$ASPL)
sd_aspl_fruit <- sd(bootstrap_fruit$ASPL)

# Same steps for stone
# ...

# Combine CC results
cc_data <- data.frame(
  CC = c(bootstrap_stone$CC, bootstrap_fruit$CC),
  Network = factor(rep(c("Stone", "Fruit"), each = 1000))
)

# Combine ASPL results
aspl_data <- data.frame(
  ASPL = c(bootstrap_stone$ASPL, bootstrap_fruit$ASPL),
  Network = factor(rep(c("Stone", "Fruit"), each = 1000))
)

# Levene's Test for CC
levene_test_cc <- leveneTest(CC ~ Network, data = cc_data)

# Levene's Test for ASPL
levene_test_aspl <- leveneTest(ASPL ~ Network, data = aspl_data)

# Welch’s t-test for ASPL
welch_test_aspl <- t.test(bootstrap_stone$ASPL, bootstrap_fruit$ASPL, var.equal = FALSE)

# Welch’s t-test for CC
welch_test_cc <- t.test(bootstrap_stone$CC, bootstrap_fruit$CC, var.equal = FALSE)


### CHI-SQUARED TESTS

# Ensure variables are treated as factors
fruit_map$Category <- as.factor(fruit_map$Category)
fruit_map$Cluster <- as.factor(fruit_map$Cluster)

# Create a contingency table
fruit_contingency_table <- table(fruit_map$Category, fruit_map$Cluster)

# Perform chi-square test
fruit_chi_test_result <- chisq.test(fruit_contingency_table, simulate.p.value = TRUE, B = 2000)

# Calculate Cramér's V using vcd
fruit_cramers_v <- assocstats(fruit_contingency_table)$cramer

# Same steps for stone
# ...
