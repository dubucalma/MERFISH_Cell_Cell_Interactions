### Cell Types proportions per layer:

library(ggplot2)
library(gridExtra)
library(dplyr)

## Load cell_infos for the different samples:
samples = c("pSTG1", "pSTG2", "pSTG3", "pSTG4", "pSTG5")
#samples = c("aSTG1", "aSTG2", "aSTG3")

## FILES WITH TRUE LAYERS OBTAINED WITH K-NEIGHBORS:
cell_infos_list <- list()
for (i in samples){
  cell_infos <- read.csv(paste0("path/to/cell_information/Layers_annotated/", i, "_annotated_cell_infos.tsv"), sep = ",")
  cell_infos_list[[i]] <- cell_infos
}


## Cell type content wihtin each layer: ONE PLOT PER SAMPLE
list_graph_distribution_per_layer<- list()
for (i in samples){
  
  cell_infos <- cell_infos_list[[i]]
  
  cell_infos$layer <- factor(cell_infos$layer)
  
  
  # Compute the proportions of each predicted.id for each bin
  prop_table <- lapply(levels(cell_infos$layer), function(lay) {
    subset <- cell_infos[cell_infos$layer == lay, ]
    prop <- prop.table(table(subset$predicted.id))
    prop_df <- data.frame(predicted.id = names(prop), percentage = prop, layer = lay)
    return(prop_df)
  })
  
  # Combine in one Data Frame
  prop_df <- do.call(rbind, prop_table)
  
  names(prop_df) <- c("predicted.id", "percentage", "layer")
  prop_df$percentage <- NULL
  names(prop_df) <- c("predicted.id", "percentage", "layer")
  
  ordered_layers <- unique(prop_df$layer)
  ordered_layers <- ordered_layers[order(as.numeric(ordered_layers))]
  
  prop_df$layer <- factor(prop_df$layer, levels = ordered_layers)
  
  # Plot the barplots
  p <- ggplot(prop_df, aes(x = layer, y = percentage, fill = predicted.id)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Percentage of each cell type within each layer", i),
         x = "Layer", y = "Percentage") +
    scale_fill_discrete(name = "Predicted ID") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  
  list_graph_distribution_per_layer[[i]] <- p
  
}


## TO HAVE THE EXACT SAME COLORS EVERY TIME: Colormap creation

all_predicted_ids <- sort(unique(unlist(lapply(cell_infos_list, function(cell_infos) {
  unique(cell_infos$predicted.id)
}))))

num_colors <- length(all_predicted_ids)
color_palette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set3"))(num_colors)  # Générer assez de couleurs

names(color_palette) <- all_predicted_ids


## ONE PLOT PER SAMPLE: 

for (i in samples) {
  
  cell_infos <- cell_infos_list[[i]]
  
  cell_infos$layer <- factor(cell_infos$layer)
  
  # Compute the proportions of each predicted.id for each bin
  prop_table <- lapply(levels(cell_infos$layer), function(lay) {
    subset <- cell_infos[cell_infos$layer == lay, ]
    prop <- prop.table(table(subset$predicted.id))
    prop_df <- data.frame(predicted.id = names(prop), percentage = prop, layer = lay)
    return(prop_df)
  })
  
  # # Combine in one Data Frame
  prop_df <- do.call(rbind, prop_table)
  
  # Add missing predicted_ids with a proportion of 0
  prop_df <- merge(expand.grid(predicted.id = all_predicted_ids, layer = unique(prop_df$layer)), 
                   prop_df, all.x = TRUE)
  prop_df$percentage.Freq[is.na(prop_df$percentage.Freq)] <- 0  
  
  ordered_layers <- unique(prop_df$layer)
  ordered_layers <- ordered_layers[order(as.numeric(ordered_layers))]
  
  prop_df$layer <- factor(prop_df$layer, levels = ordered_layers)
  
  prop_df$predicted.id <- factor(prop_df$predicted.id, levels = all_predicted_ids)
  
  # Plot the barplot with a fixed color panel
  p <- ggplot(prop_df, aes(x = layer, y = percentage.Freq, fill = predicted.id)) +
    geom_bar(stat = "identity", position = "stack") +  
    labs(title = paste("Percentage of each cell type within each layer", i),
         x = "Layer", y = "Percentage") +
    scale_fill_manual(values = color_palette) +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = "white"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  
  list_graph_distribution_per_layer[[i]] <- p
}


grid.arrange(grobs = list_graph_distribution_per_layer, ncol = 3)



#### AGGREGATION FOR ALL SAMPLES IN ONE PLOT per pSTG or aSTG:
combined_cell_infos <- data.frame()

# Combine all samples
for (i in samples) {
  cell_infos <- cell_infos_list[[i]]
  cell_infos$sample <- i  # Ajouter une colonne pour identifier le sample
  combined_cell_infos <- rbind(combined_cell_infos, cell_infos)
}

# Aggregate the cells per cell types, per sample:
sum_table <- aggregate(. ~ layer + predicted.id, data = combined_cell_infos, FUN = length)

colnames(sum_table)[colnames(sum_table) == "sample"] <- "count"

# Proportions per layer
sum_table <- sum_table %>%
  group_by(layer) %>%
  mutate(percentage = count / sum(count)) %>%
  ungroup()

# Ordonner les layers numériquement si ce n'est pas déjà le cas
ordered_layers <- unique(sum_table$layer)
ordered_layers <- ordered_layers[order(as.numeric(ordered_layers))]

sum_table$layer <- factor(sum_table$layer, levels = ordered_layers)

# Save proportions table as CSV --> SEE LOCAL CODE

# Barplot
p <- ggplot(sum_table, aes(x = layer, y = percentage, fill = predicted.id)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Proportion of each cell type across pSTG samples by layer",
       x = "Layer", y = "Proportion") +
  scale_fill_discrete(name = "Predicted ID") +
  scale_fill_manual(values = color_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = paste0("path/to/cell_information/Plots/Human_pSTG_cell_content_sept2025.pdf"), 
       plot = p, device = "pdf", width = 15, height = 10)


  