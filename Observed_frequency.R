## Computation of OBSERVED Frequency: 
library(dplyr)
library(dbscan)

## DATA NEEDED: 
#study_dir = ""

samples = c("pSTG1", "pSTG2", "pSTG3", "pSTG4", "pSTG5")
eps = 15 

## FUNCTION: 
calculate_contact_frequency <- function(cell_positions, cell_types, radius) {

  '''Given cell coordinates and their predicted cell types,
    this function computes neighbor contact frequencies within
    a fixed radius (eps = 15 µm).'''

  nn <- frNN(x = cell_positions[, c("x_final", "y_final")], eps = eps)
  # Create a list of tables for each cell
  x_cell_types <- purrr::map(nn$id, ~ as.factor(cell_types$predicted.id)[.x] %>% table())
  # Combine the tables into a matrix
  nn_matrix <- do.call(rbind, x_cell_types)
  nn_df <- as.data.frame(nn_matrix)
  df <- tibble(CellType = cell_types$predicted.id, as.data.frame(nn_matrix))
  result_df <- df %>%
    group_by(CellType) %>%
    summarise(across(everything(), sum))
  
  return(result_df)
}

tibble_data_frame_format_matrix <- function(tibble){

  '''Convert a tibble of contact frequencies (with a `CellType` column)
      into a clean numeric matrix.'''

  contact_data_df <- as.data.frame(tibble)
  row.names(contact_data_df) <- tibble$CellType
  contact_data_df <- contact_data_df[, -which(colnames(contact_data_df) == "CellType")]
  # Sort row names alphabetically
  sorted_row_names <- sort(rownames(contact_data_df))
  contact_data_df <- contact_data_df[sorted_row_names, ]
  # Sort column names alphabetically
  sorted_col_names <- sort(colnames(contact_data_df))
  contact_data_df <- contact_data_df[, sorted_col_names]
  
  data_matrix <- as.matrix(contact_data_df)
  return(data_matrix)
}

for (i in samples){

  cell_infos <- read.table(paste0("path/to/cell_information/", i, "superclass_layer_annotated_cell_infos.tsv"), sep = "\t")
  
  list_bins <- list()
  
  list_layers = c("upper", "deeper")
  
  for (bin in list_layers){
    print(bin)
    subset_bin <- subset(cell_infos, slided_window == bin)
    subset_bin$x_final <- as.numeric(subset_bin$x_final)
    subset_bin$y_final <- as.numeric(subset_bin$y_final)
    observed_frequency <- calculate_contact_frequency(subset_bin, subset_bin, eps)
    observed_frequency_matrix <- tibble_data_frame_format_matrix(observed_frequency)
    list_bins[[bin]] <- observed_frequency_matrix

  }
  
  saveRDS(list_bins, paste0("path/to/cell_information/Observed_frequencies/", i, "_layers_cells_observed_freq.rds"))
  
}


