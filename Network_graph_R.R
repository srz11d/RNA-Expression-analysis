
# Load required library
library(igraph)

# Read the reshaped data
reshaped_data <- read.csv("Shaped_network.csv")

# Create a graph object
graph <- graph_from_data_frame(reshaped_data[, c("ALT", "Genes", "Type", "Spearman")])

# Plot the graph
plot(graph, edge.width = abs(E(graph)$Spearman) * 3, edge.color = ifelse(E(graph)$Spearman > 0, "blue", "red"), vertex.label.cex = 0.8, vertex.label.dist = 2)


###
# Read the reshaped data
reshaped_data <- read.csv("Shaped_network.csv")

# Create a graph object
graph <- graph_from_data_frame(reshaped_data[, c("ALT", "Genes", "Type", "Spearman")])

# Set up node attributes (e.g., shape and color)
V(graph)$type <- ifelse(V(graph)$name == "ALT", "Type", "Spearman")
V(graph)$color <- ifelse(V(graph)$type == "ALT", "red", "lightblue")

# Set up edge attributes (e.g., color)
E(graph)$edge_color <- ifelse(E(graph)$Spearman > 0, "magenta", "gray")

# Customize layout
layout <- layout_with_fr(graph)

# Plot the graph with customizations
plot(
  graph,
  layout = layout,
  vertex.size = 8,
  vertex.label.cex = 0.6,
  vertex.color = V(graph)$color,
  vertex.shape = ifelse(V(graph)$type == "ALT", "rectangle", "circle"),
  edge.width = abs(E(graph)$Spearman) * 3,
  edge.color = E(graph)$edge_color,
  main = "Network Graph of ALT Cells and Other Cell Types"
)