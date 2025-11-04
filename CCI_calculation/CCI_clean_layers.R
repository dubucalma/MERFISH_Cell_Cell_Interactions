## CCI for layers - bins:


#############################
## FUNCTIONS TO FORMAT data: 
############################

## Format data + Just take the upper part of the data.frame: 
tibble_data_frame_format_matrix <- function(tibble){
  contact_data_df <- as.data.frame(tibble)
  row.names(contact_data_df) <- tibble$CellType
  contact_data_df <- contact_data_df[, -which(colnames(contact_data_df) == "CellType")]
  # Sort row names alphabetically
  sorted_row_names <- sort(rownames(contact_data_df))
  contact_data_df <- contact_data_df[sorted_row_names, ]
  # Sort column names alphabetically
  sorted_col_names <- sort(colnames(contact_data_df))
  contact_data_df <- contact_data_df[, sorted_col_names]
  
  # Convert the data frame to a matrix
  data_matrix <- as.matrix(contact_data_df)
  return(data_matrix)
}

## For a given sample, creation of a list of symmetric matrices:  
Contact_matrix_list_creation <- function(sample, bin){
  tibble_list <- get(paste0(sample, "_randomized_contacts_layers"))[[bin]]
  Contacts_matrix_list <- list()
  for (i in seq_along(tibble_list)){
    matrix <- tibble_data_frame_format_matrix(tibble_list[[i]])
    Contacts_matrix_list[[i]] <- matrix
  }
  
  return (Contacts_matrix_list)
}

## For a given sample, creation of a list of symmetric matrices:  
Contact_matrix_list_creation <- function(sample){
  tibble_list <- get(paste0(sample, "_randomized_contacts_layers"))
  Contacts_matrix_list <- list()
  for (i in seq_along(tibble_list)){
    matrix <- tibble_data_frame_format_matrix(tibble_list[[i]])
    Contacts_matrix_list[[i]] <- matrix
  }
  
  return (Contacts_matrix_list)
}

#Transform the list of tibble in a list of symmetric matrices: 
for (i in samples){
  name <- paste0(i, "_randomized_contacts_matrix_list") 
  
  list_permut <- get(paste0(i, "_randomized_contacts_layers")) 
  
  bins_matrix_list <- list()
  
  layer_superclasses = c("deeper", "upper")
  
  for (bin in layer_superclasses){
    bin_list_permut <- list_permut[[bin]]
    matrix_list <- Contact_matrix_list_creation(i, bin)
    bins_matrix_list[[bin]] <- matrix_list
  }
  
  assign(name, bins_matrix_list)
}


############################
## FUNCTIONS for Statistics
############################

# Function for each coordinate, through all matrices:
calculate_mean_matrice <- function(list_of_matrices) {
  # Initialize a matrix to store the mean values for each coordinate
  n <- nrow(list_of_matrices[[1]])
  m <- ncol(list_of_matrices[[1]])
  mean_values <- matrix(0, nrow = n, ncol = m)
  
  # Calculate the mean value for each coordinate across all matrices
  for (i in 1:length(list_of_matrices)) {
    mean_values <- mean_values + list_of_matrices[[i]]
  }
  mean_values <- mean_values / length(list_of_matrices)
  
  return(mean_values)
} # Expected frequencies, null mean

calculate_sd_matrice <- function(list_of_matrices, mean_values) {
  # Initialize a matrix to store the standard deviation values for each coordinate
  n <- nrow(list_of_matrices[[1]])
  m <- ncol(list_of_matrices[[1]])
  sd_values <- matrix(0, nrow = n, ncol = m)
  
  # Step 2: Calculate the sum of squared differences
  for (i in 1:length(list_of_matrices)) {
    sd_values <- sd_values + (list_of_matrices[[i]] - mean_values)^2
  }
  
  # Step 3: Divide the sum by (n - 1)
  sd_values <- sqrt(sd_values / (length(list_of_matrices) - 1))
  
  return(sd_values)
} # null sd

calculate_fold_change_matrice <- function(measured_mat, expected_mat) {
  fold_change <- (measured_mat - expected_mat) / expected_mat
  return(fold_change)
}

calculate_z_score_matrice <- function(observed_mat, null_mean_mat, null_sd_mat) {
  z_score <- (observed_mat - null_mean_mat) / null_sd_mat
  return(z_score)
}

calculate_p_value_matrice <- function(zscore_matrix) {
  # Apply the pnorm function to each element of the matrix
  p_values <- apply(zscore_matrix, c(1, 2), function(x) 2 * (1 - pnorm(abs(x))))
  
  return(p_values)
}

calculate_FDR_matrix <- function(p_values_matrix) {
  """ This function aims at filtering the matrix with a threshold on the corrected Pvalue only """

  n <- nrow(p_values_matrix)
  
  # Extract upper triangle (including diagonal) of the p-values matrix
  upper_triangle <- p_values_matrix
  upper_triangle[lower.tri(upper_triangle)] <- "NA"
  
  # Compute the FDR correction for upper triangle using p.adjust
  fdr_adjusted_upper <- p.adjust(upper_triangle[upper.tri(upper_triangle)], method = "fdr")
  
  # Reconstruct the symmetric matrix with FDR-adjusted p-values
  fdr_matrix <- p_values_matrix
  fdr_matrix[upper.tri(fdr_matrix)] <- fdr_adjusted_upper
  
  # Fill the lower triangle with (mirrored) NA values
  fdr_matrix[lower.tri(fdr_matrix)] <- NA
  
  return(fdr_matrix)
}

calculate_FDR_BH_matrix <- function(p_values_matrix) {
  """ This function aims at filtering the matrix with a thershold on the BH - corrected Pvalue  """

  n <- nrow(p_values_matrix)
  
  # Extract upper triangle (including diagonal) of the p-values matrix
  upper_triangle <- p_values_matrix
  upper_triangle[lower.tri(upper_triangle)] <- NA
  
  # Compute the BH-adjusted p-values for the upper triangle
  bh_adjusted_upper <- p.adjust(upper_triangle[upper.tri(upper_triangle)], method = "BH")
  
  # Reconstruct the symmetric matrix with BH-adjusted p-values
  fdr_matrix_bh <- p_values_matrix
  fdr_matrix_bh[upper.tri(fdr_matrix_bh)] <- bh_adjusted_upper
  
  # Fill the lower triangle with mirrored values
  fdr_matrix_bh[lower.tri(fdr_matrix_bh)] <- fdr_matrix_bh[upper.tri(fdr_matrix_bh)]
  
  return(fdr_matrix_bh)
} ## Correction with Benjamini and Hochberg procedure.

filter_fdr_significant <- function(fdr_matrix){
  fdr_matrix_filtered <- fdr_matrix
  fdr_matrix_filtered[fdr_matrix >= 0.01] <- NA
  return(fdr_matrix_filtered)
}

extract_and_recreate_matrices <- function(fold_change_matrix, fdr_matrix, fold_change_threshold, fdr_threshold) {
  significant_indices <- which(fold_change_matrix > fold_change_threshold & fdr_matrix < fdr_threshold, arr.ind = TRUE)

  filtered_fold_change <- fold_change_matrix[significant_indices]
  filtered_fdr <- fdr_matrix[significant_indices]
  
  # Create an empty matrix 
  recreated_fold_change <- matrix(NA, nrow = nrow(fold_change_matrix), ncol = ncol(fold_change_matrix))
  recreated_fdr <- matrix(NA, nrow = nrow(fdr_matrix), ncol = ncol(fdr_matrix))
  
  # Fill the matrix with significant data
  recreated_fold_change[significant_indices] <- filtered_fold_change
  recreated_fdr[significant_indices] <- filtered_fdr
  
  # Name Rows / columns
  rownames(recreated_fold_change) <- rownames(fold_change_matrix)
  colnames(recreated_fold_change) <- colnames(fold_change_matrix)
  rownames(recreated_fdr) <- rownames(fdr_matrix)
  colnames(recreated_fdr) <- colnames(fdr_matrix)
  
  return(list(recreated_fold_change = recreated_fold_change, recreated_fdr = recreated_fdr))
}

count_non_na <- function(mat){
  count <- sum(!is.na(mat))
  return(count)
}



#################
#### PARAMETERS :
##################
samples = c("pSTG1", "pSTG2", "pSTG3", "pSTG4", "pSTG5")


# Loading 1000 permutation results for all sample: 
for (i in samples){
  name <- paste0(i, "_randomized_contacts_layers") 
  permut <- readRDS(paste0("path/to/cell_information/CCI_randomization/", i, "_1000_Permuted_mean_frequencies.rds"))
  assign(name, permut)
}


#################
# FORMAT DATA
#################

###### FOR WHOLE SLICE (NO BIN, NO LAYER)
for (i in samples){
  name <- paste0(i, "_randomized_contacts_matrix_list")
  tibble_list <- get(paste0(i, "_randomized_contacts_layers"))
  matrix_list <- Contact_matrix_list_creation(i)
  assign(name, matrix_list)
}


####################################### 
## CCI results for each sample PER BIN: 
######################################

for (i in samples){

  observed_frequency_total <- readRDS(paste0("path/to/cell_information/Observed_frequencies/", i, "_layers_cells_observed_freq.rds"))
  
  #nb_bins <- length(observed_frequency_total) # FOR LAYERS, supposed to have just 2, "Upper and deeper"
  
  Results_all_sample <- list()
  
  
  layer_superclasses = c("deeper", "upper")

  
  for (bin in layer_superclasses){
    observed_frequency <- observed_frequency_total[[bin]]
       
    expected_frequency <- calculate_mean_matrice(get(paste0(i, "_randomized_contacts_matrix_list"))[[bin]])
    sd <- calculate_sd_matrice((get(paste0(i, "_randomized_contacts_matrix_list")))[[bin]], expected_frequency)
    
    intersection = intersect(colnames(observed_frequency), colnames(expected_frequency))
    expected_frequency_clean = expected_frequency[intersection, intersection]
    observed_frequency_clean = observed_frequency[intersection, intersection]
    sd_clean <- sd[intersection, intersection]
    
    fold_change <- calculate_fold_change_matrice(observed_frequency_clean, expected_frequency_clean)
    zscore <- calculate_z_score_matrice(observed_frequency_clean, expected_frequency_clean, sd_clean)
    pvalue <- calculate_p_value_matrice(zscore)
    fdr <- calculate_FDR_matrix(pvalue)
    fdr_significant <- filter_fdr_significant(fdr)
    
    # Filter on significant values.
    filtered_results <- extract_and_recreate_matrices(fold_change, fdr, fold_change_threshold = 1.5, fdr_threshold = 0.05)
    
    count_interac_strong <- count_non_na(filtered_results$recreated_fdr) # To count strong interactions...
    
    
    results_list <- list(
      expected_frequency = expected_frequency,
      standard_deviation = sd,
      fold_change = fold_change,
      zscore = zscore,
      pvalue = pvalue,
      fdr = fdr,
      fdr_significant = fdr_significant,
      recreated_fold_change = filtered_results$recreated_fold_change, ## Matrix with only both sides significant contacts
      recreated_fdr = filtered_results$recreated_fdr, ## Same, fdr values
      counts_interac <- count_interac_strong ## to be called as [[10]].
    )
    
    Results_all_sample[[bin]] <- results_list
    write.table(fold_change, paste0("path/to/cell_information/Fold_changes/", i, "/fold_change_", bin, ".tsv"), sep = "\t")
    
  }
  
  name <- paste0(i, "_CCI_recap_list")
  assign(name, Results_all_sample)
  saveRDS(Results_all_sample, paste0("path/to/cell_information/Summary_CCI/", i, "_CCI_per_layer_recap.rds"))

}