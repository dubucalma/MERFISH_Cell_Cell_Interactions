
## Interactions analysis (cell type specific):

###########
##LIBRARIES
###########
library(dbscan)
library(dplyr)
library(future.apply)
plan(multisession, workers=20)

##############
###PARAMETERS 
#############
args = commandArgs(trailingOnly = TRUE)

samples = c("pSTG1", "pSTG2", "pSTG3", "pSTG4", "pSTG5")
sample = samples[as.integer(args)]

#study_dir = 

radius = 100
eps = 15
perm_count = 1000

#################
#.   FUNCTIONS : 
#################

Cell_Neighbors_function <- function(sample){

    """ Given a sample ID, looks up `cell_infos_<sample>` (must contain `predicted.id`)
  #' and `coordinates_<sample>` (must contain columns `cell`, `x`, `y`) in the
  #' environment, finds neighbors within `eps` micrometers using `dbscan::frNN`,
  #' and aggregates counts by focal cell type (rows) × neighbor cell type (cols)."""
  
  # Parameters for the cell: 
  cell_infos <- get(paste0("cell_infos_", sample))
  coordinates <- get(paste0("coordinates_", sample))
  rownames(coordinates) <- coordinates$cell
  coordinates <- select(coordinates, - "cell")
  coordinates <- as.matrix(coordinates)
  
  
  ## Finding the closest neighbors: 15 micrometers
  eps <- 15
  nn <- frNN(x= coordinates, eps = eps)
  
  ## Neighborhood Count Matrix: 
  x_cell_types<- purrr::map(nn$id, ~as.factor(cell_infos$predicted.id)[.x] %>% table())
  nn_matrix_asd<- do.call(rbind,x_cell_types)
  head(nn_matrix_asd) # Columns = neighbor (contact) cell types from `predicted.id`; rows = focal cells
  # For each focal cell i, entry (i, t) = number of contacts (neighbors within eps) with cell type t
  
  #Computation of the contact or proximity frequency between two cell-types A/B: 
  predictions_simple = select(seurat_subset@meta.data, "predicted.id")
  nn_df_asd <- as.data.frame(nn_matrix_asd)

  # Gather all data in a tibble
  df_asd <- tibble(CellType = predictions_simple, as.data.frame(nn_matrix_asd))
  
  result_df_asd <- df_asd %>%
    group_by(CellType) %>%
    summarise(across(everything(), sum))
  print(result_df_asd)
} 

# Function to randomize cell positions within a given radius
randomize_positions <- function(x, y, radius) {
  theta <- runif(length(x), 0, 2 * pi)
  r <- runif(length(x), 0, radius)
  x_rand <- x + r * cos(theta)
  y_rand <- y + r * sin(theta)
  return(data.frame(x = x_rand, y = y_rand))
}

calculate_contact_frequency <- function(cell_positions, cell_types, radius) {

  '''Given cell coordinates and their predicted cell types,
    this function computes neighbor contact frequencies within
    a fixed radius (eps = 15 µm).'''
  
  nn <- frNN(x = cell_positions[, c("x", "y")], eps = eps)  
  # Create a list of tables for each cell
  x_cell_types <- purrr::map(nn$id, ~ as.factor(cell_types$predicted.id)[.x] %>% table())
  nn_matrix <- do.call(rbind, x_cell_types)
  nn_df <- as.data.frame(nn_matrix)
  df <- tibble(CellType = cell_types$predicted.id, as.data.frame(nn_matrix))
  result_df <- df %>%
    group_by(CellType) %>%
    summarise(across(everything(), sum))
  
  return(result_df)
} 

Apply_freq_to_permut <- function(i){

  '''Perform one permutation (i):
      - calls randomize_positions() to jitter cell coordinates
      - calls calculate_contact_frequency() to compute the contact matrix
      - prints runtime and returns the result '''

  start_time <- Sys.time()
  coordinates_rand <- randomize_positions(as.numeric(as.data.frame(cell_infos)$x_final), as.numeric(as.data.frame(cell_infos)$y_final), radius)
  
  merged_data_rand <- data.frame(
    x = coordinates_rand$x,
    y = coordinates_rand$y,
    predicted_id = cell_infos$predicted.id[rownames(coordinates_rand)]
  )
  
  perm_results <- calculate_contact_frequency(merged_data_rand, cell_infos, eps)
  
  end_time <- Sys.time()  
  time_diff <- difftime(start_time, end_time, units = "secs") 
  print(paste("Iteration", i, "took", round(time_diff,3) , "seconds"))
  return(perm_results)
}


## Function with permutation f randomization + calulate mean every time: 
Compute_frequency_randomization <- function(sample, perm_count) {

  '''Run multiple permutations for a given sample:
        - calls Apply_freq_to_permut() perm_count times (in parallel with future_lapply)
        - collects and returns the list of all contact frequency matrices '''
  
  cell_infos <- get(paste0("cell_infos_", sample))
  
  all_perm_results <- vector("list", perm_count)
  all_perm_results <- future_lapply(seq(1,perm_count,1), Apply_freq_to_permut, future.seed = TRUE)
  
  return(all_perm_results)
}


name <- paste0("cell_infos_", sample)
cell_infos <- read.table(paste0("path/to/cell_information/", sample, "_cell_infos.tsv"), sep = ",")
assign(name, cell_infos)

name_bis <- paste0("perm_results_", sample)
perm_results <- Compute_frequency_randomization(sample, perm_count)
  
saveRDS(perm_results, paste0("path/to/cell_information/CCI_randomization/", i, "_1000_Permuted_mean_frequencies.rds"))
assign(name_bis, perm_results)



