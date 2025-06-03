# Title: Testing the effect of scale on species-landscape relationships
# Author: Herbert Leavitt 
# Date: February 17, 2025

# Install and load required packages
# Install 'pscl' package if not already installed
if (!require(pscl)) install.packages("pscl")
# Load necessary libraries
library(statmod)    # For statistical modeling functions
library(scales)     # For percentage formats in ggplots
library(tidyverse)  # For data manipulation and visualization
library(spdep)      # For spatial dependence and spatial regression
library(gam)        # For generalized additive models
library(mgcv)       # For advanced GAMs with automatic smoothing parameter selection
library(ggplot2)    # For data visualization
library(stringr)
formula_table <- list() # Initiate formula table, don't run this line again. 


run_scale = "smallscale" # should be "smallscale" for drone imagery of "satscale" for satellite" 

outlier_table <- data.frame(species = c("PALSP"), site_date_key = c("05244dc3"), stringsAsFactors = FALSE) #Add points designated as outliers. 

######## Function Definitions ###########
# Back-transform function for plotting
back_transform <- function(value, var, stats) {
  mean_val <- stats[stats$variable == var, "mean"]
  sd_val <- stats[stats$variable == var, "sd"]
  (value * sd_val) + mean_val
}

# Function to check if a model contains any highly correlated pairs of predictors
contains_correlated_pair <- function(predictors, correlated_pairs) {
  # Iterate over each pair of correlated predictors
  for (i in 1:nrow(correlated_pairs)) {
    # Check if both predictors in the pair are included in the current model
    if (all(correlated_pairs[i, ] %in% predictors)) {
      return(TRUE)  # Model contains a highly correlated pair
    }
  }
  return(FALSE)  # No highly correlated pairs in the model
}

# Function to create a table with one variable sequenced and others held at their mean
create_sequence_table <- function(data, sequence_var, length_out = 100) {
  # Check if the sequence_var exists in data and is numeric
  if (!sequence_var %in% names(data)) {
    stop("The specified sequence_var does not exist in the data.")
  }
  
  if (!is.numeric(data[[sequence_var]])) {
    stop("The specified sequence_var is not numeric.")
  }
  
  # Calculate means for numeric columns only and create sequence
  data %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE), .names = "{col}")) %>%
    slice(rep(1, length_out)) %>%
    mutate(!!sequence_var := seq(min(data[[sequence_var]], na.rm = TRUE),
                                 max(data[[sequence_var]], na.rm = TRUE),
                                 length.out = length_out))
}


# Main function for model compilation and selection. The 'model selection' function is nested within this function. The purpose of this function is to run the preliminary 
# data wrangling for each iteration of buffer and edge habitat. This includes checking for correlated variables and standardizing all variables. After 'model_selection' is run, it 
# tidies up the outputs into a usable table that's retured. 
# Arguments:
# - comtab: Dataframe community table abundance data
# - species_code: Code for the species being modeled
# - buffers: Numeric vector of buffer distances to iterate over
# - edges: Numeric vector of edge distances to iterate over
# - aic_threshold: AIC threshold for model selection (default is 2)
# - correlation_threshold: Threshold for defining highly correlated pairs (default is 0.8)
# - scale: run scale that was selected earlier 
model_compile <- function(comtab, species_code, buffers, edges, aic_threshold = 2, correlation_threshold = 0.8, scale = run_scale) {
  # Initialize the data frame to store results
  r.tab <- data.frame(
    edge = numeric(),
    buffer = numeric(),
    sp = character(),
    best_model = character(),
    aic = numeric(),
    r.sq = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop over each combination of buffer size and edge distance
  for (buf in buffers) {
    for (edge in edges) {
      
      if( scale == "satscale"){
        # Construct the filename for habitat data
        habx <- paste("google2022_", "edge", edge, "_buf", buf, sep = "")
        
        # Load habitat data
        hab <- read.csv(file.path("./landscape_analysis", "output", "satscale", paste(habx, ".csv", sep = "")))
      }else{
        habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
        
        # Load habitat data
        hab <- read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx, ".csv", sep = "")))
      }
      # Calculate percentage of edge habitat
      hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
      
      # Merge habitat data with community composition data
      pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key')
      
      # Generate correlation matrix and identify correlated pairs
      numeric_data <- pf.env[, c("edge_perc", "edge_l.mangrove", "land_water_ratio")]
      all_predictors <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
      all_predictors_and_interactions <- c("edge_perc", "edge_l.mangrove", "land_water_ratio", "man_lwr", "man_edge")
      correlation_matrix <- cor(numeric_data, use = "complete.obs")
      standardization_stats <- data.frame(variable = all_predictors, mean = numeric(length(all_predictors)), sd = numeric(length(all_predictors)))
      
      for (col in all_predictors) {
        standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
        standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
        pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
      }
      # Identify highly correlated pairs (absolute correlation > threshold)
      high_corr <- which(abs(correlation_matrix) > correlation_threshold & upper.tri(correlation_matrix), arr.ind = TRUE)
      correlated_pairs <- data.frame(
        predictor1 = rownames(correlation_matrix)[high_corr[, 1]],
        predictor2 = colnames(correlation_matrix)[high_corr[, 2]]
      )
      
      # Run model selection for the current combination
      selected_models <- model_selection(pf.env, response_var = species_code, all_predictors = all_predictors_and_interactions, correlated_pairs, aic_threshold, buf, edge, scale)
      
      # Collect formula and AIC for each model within the AIC threshold
      for (i in seq_along(selected_models)) {
        model_info <- list(
          formula = formula(selected_models[[i]]),
          aic = AIC(selected_models[[i]]),
          r.sq = summary(selected_models[[i]])$r.sq
        )
        
        # Create a new row with buffer, edge, species, model formula, and AIC
        newline <- data.frame(
          edge = edge,
          buffer = buf,
          sp = species_code,
          best_model = paste(deparse(model_info$formula), collapse = " "),
          aic = model_info$aic,
          r.sq =  model_info$r.sq,
          stringsAsFactors = FALSE
        )
        
        # Append the new row to the results table
        r.tab <- rbind(r.tab, newline)
      }
    }
  }
  
  return(r.tab) # retuen results table
}

# Model selection function. This is nested within the 'model_compile' function. It takes all the data from the iteration of buffer and edge determined from 'model_compile' along with information on
# possible correlated pairs. Then it selects models within +/- 2 delta AIC for each species.
# Arguments:
# - data: Data frame containing predictors and response variables
# - response_var: Variable name for the response variable
# - all_predictors: Vector of predictor variable names
# - correlated_pairs: Data frame of correlated predictor pairs
# - aic_threshold: AIC threshold for selecting models (default is 2)
model_selection <- function(data, response_var, all_predictors, correlated_pairs, aic_threshold = 2, buf, edge, scale = scale) {
  best_aic <- Inf
  models_all <- list()  # Store models within the AIC threshold
  all_models_info <- data.frame(formula = character(), aic = numeric(), stringsAsFactors = FALSE)  # Store all tested models
  
  # Iterate over all subsets of predictors to get every possible combination
  for (subset_size in 1:length(all_predictors)) {
    predictor_subsets <- combn(all_predictors, subset_size, simplify = FALSE)
    valid_subsets <- Filter(function(predictor_subsets) {
      subset <- unlist(predictor_subsets)
      if ("man_edge" %in% subset & !("edge_l.mangrove" %in% subset & "edge_perc" %in% subset)) {
        return(FALSE)  # Remove this subset
      }
      if ("man_lwr" %in% subset & !("edge_l.mangrove" %in% subset & "land_water_ratio" %in% subset)) {
        return(FALSE)  # Remove this subset
      }
      return(TRUE)  # Keep subset if it meets the conditions
      
    }, predictor_subsets)
    
    #for loop now runs for every possible combination of valid model configurations
    for (subset_predictors in valid_subsets) {
      # Skip subsets containing highly correlated pairs
      if (contains_correlated_pair(subset_predictors, correlated_pairs)) next
      
      # Create formula terms with smooth terms for each predictor
      # Add interaction terms if applicable
      if ("man_lwr" %in% subset_predictors & "man_edge" %in% subset_predictors) {
        subset_predictors <- setdiff(subset_predictors, c("man_lwr", "man_edge"))
        formula_terms <- paste(sprintf("s(%s, k = 4)", subset_predictors), collapse = " + ")
        formula_terms <- paste(formula_terms, " + ti(edge_l.mangrove, land_water_ratio, k = 3)+ ti(edge_l.mangrove, edge_perc, k = 3)")
      } else if ("man_lwr" %in% subset_predictors) {
        subset_predictors <- setdiff(subset_predictors, "man_lwr")
        formula_terms <- paste(sprintf("s(%s, k = 4)", subset_predictors), collapse = " + ")
        formula_terms <- paste(formula_terms, " + ti(edge_l.mangrove, land_water_ratio, k = 3)")
      } else if ("man_edge" %in% subset_predictors) {
        # Interaction between edge length of mangrove and edge percentage
        subset_predictors <- setdiff(subset_predictors, "man_edge")
        formula_terms <- paste(sprintf("s(%s, k = 4)", subset_predictors), collapse = " + ")
        formula_terms <- paste(formula_terms, " + ti(edge_l.mangrove, edge_perc, k = 3)") 
      }else{
        formula_terms <- paste(sprintf("s(%s, k = 4)", subset_predictors), collapse = " + ")
        
      }

      # Construct the full formula and fit the GAM model
      formula <- as.formula(paste(response_var, "~", formula_terms))
      model <- gam(formula, data = data, family = tw(link = "log"), method = "ML")
      model_aic <- AIC(model)
      model_rsq <- summary(model)$r.sq
      models_all <- c(models_all, list(model))  
      
      # Store the model formula and AIC
      all_models_info <- rbind(all_models_info, data.frame(formula = paste(deparse(formula), collapse = " "), 
                                                           aic = model_aic,
                                                           r.sq = model_rsq,
                                                           stringsAsFactors = FALSE))
      # all model info gets written into the repository 
      write.csv(all_models_info, file.path("species_analysis", "stat_intermediates", "species_aic", paste0(species_code,"_", scale, "edge", edge, "buf", buf, ".csv")))
      
      # Update the list of best models based on AIC
      if (model_aic < best_aic) {
        best_aic <- model_aic
      } 
    }
  }
  
  # Extract AIC values
  aic_values <- sapply(models_all, AIC)
  
  #extract models within threshold of "best" AIC 
  models_within_aic <- models_all[aic_values <= best_aic + aic_threshold]
  # Save all tested models and their AIC scores
  write.csv(all_models_info, file.path("species_analysis", "stat_intermediates", "species_aic", paste0(scale,"", ".csv")))
  
  return(models_within_aic)
}

# Now knowing the correct model for the species at given scale, this function runs the model across all scales and produces 
# Rsq plots and partial effect plots 
# Arguments 
# - comtab: community abundance table 
# - species_code: character string indicating the species to be modeled 
# - edges: vector of edge distances 
# - buffers: vector of buffer distances 
# - top_model_formula: formula of the model selected above 
# - scale: satscale or smallscale
# Function to run model across scales with dynamic handling of variables
univariate_edge <- function(comtab, species_code, edges, buffers, top_model_formula, scale = run_scale) {
  effect_data <- data.frame()
  interaction_variables <- vector()
  # Extract all terms as they appear in the formula, including interaction terms
  model_terms <- all.vars(top_model_formula)  # Gets all unique variable names
  all_formula_terms <- labels(terms(top_model_formula))  # Gets all terms as they appear in the formula
  
  # Identify interaction terms explicitly by looking for "ti(" pattern
  interaction_term <- all_formula_terms[grepl("ti\\(", all_formula_terms)]
  
  # Parse interaction terms to extract each variable involved
  cleaned_interaction_term <- gsub("ti\\(([^,]+),([^,]+).*\\)", "ti(\\1,\\2)", interaction_term)
  
  if (length(cleaned_interaction_term) != 0) {
    # Check the value of cleaned_interaction_term for specific interactions
    
    if ("ti(edge_l.mangrove, land_water_ratio)" %in% cleaned_interaction_term) {
      interaction_variables <- c(interaction_variables, "man_lwr")
    } else if (cleaned_interaction_term == "ti(edge_l.mangrove, edge_perc)") {
      interaction_variables <-  c(interaction_variables,"man_edge")
    }
  }
  
  # Initialize essential columns and dynamically add relevant columns based on the formula
  columns <- c("edge", "buffer", "sp", "r_sq", "moran_p", "cov", "dev_ex")  # Essential columns
  variable_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio",  "man_lwr", "man_edge")  # Possible predictor variables
  
  # Combine single terms and interaction terms
  
  selected_columns <- intersect(model_terms, variable_columns)
  if (length(cleaned_interaction_term) != 0) {
    selected_columns <- unique(c(selected_columns, interaction_variables))
  }
  # Add "p." prefixed columns for p-values of each selected term
  result_columns <- c(columns, paste0("p.", selected_columns))
  
  # Create the results table with dynamically selected columns
  r.tab <- as.data.frame(matrix(ncol = length(result_columns), nrow = 0))
  colnames(r.tab) <- result_columns
  
  # Loop over each combination of edge distance and buffer size - store all the data for output and plotting 
  for (edge in edges) {
    for (buf in buffers) {
      tryCatch({
        # Construct the filename for habitat data
        
        if(scale == "satscale"){
          # Construct the filename for habitat data
          habx <- paste("google2022_", "edge", edge, "_buf", buf, sep = "")
          
          # Load habitat data
          hab <- read.csv(file.path("./landscape_analysis", "output", "satscale", paste(habx, ".csv", sep = "")))
        }else{
          habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
          
          # Load habitat data
          hab <- read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx, ".csv", sep = "")))
        }
        # Calculate percentage of edge habitat
        hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
        
        # Merge habitat data with community composition data
        pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key')
        # List of all predictors
        all_predictors <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
        numeric_data <- pf.env[, c("edge_perc", "edge_l.mangrove", "land_water_ratio")]
        
        # standardize data once again
        standardization_stats <- data.frame(variable = all_predictors, mean = numeric(length(all_predictors)), sd = numeric(length(all_predictors)))
        
        for (col in all_predictors) {
          standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
          standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
          pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
        }
        
        
        
        # Fit the top model
        gam_model <-  gam(as.formula(top_model_formula), data = pf.env, family = tw(link = "log"), method = "ML")
        
        # Get GAM model summary and extract model outputs
        gam_summary <- summary(gam_model)
        r_sq <- gam_summary$r.sq
        dev_ex <- gam_summary$dev.expl
        # Set the threshold for high correlation
        correlation_threshold <- 0.8
        # Calculate residuals and plot diagnostics
        par(mfrow = c(2, 2))
        
        fitted_values <- fitted(gam_model)
        residuals_gam <- residuals(gam_model, type = "pearson")
        
        valid_indices <- is.finite(fitted_values) & is.finite(residuals_gam)
        if (any(!valid_indices)) {
          warning("Non-finite fitted values or residuals detected. Plotting will exclude these.")
          fitted_values <- fitted_values[valid_indices]
          residuals_gam <- residuals_gam[valid_indices]
        }
        
        # Plot only valid values
        par(mfrow = c(2, 2))
        
        # Q-Q Plot of Residuals
        qqnorm(residuals_gam, main = "Q-Q Plot of GAM Residuals")
        qqline(residuals_gam, col = "red", lwd = 2)
        
        # Residuals vs. Fitted Values
        plot(fitted_values, residuals_gam,
             xlab = "Fitted Values", ylab = "Residuals",
             main = "Residuals vs Fitted Values")
        abline(h = 0, col = "red", lwd = 2)
        
        # Scale-Location Plot
        plot(fitted_values, sqrt(abs(residuals_gam)),
             xlab = "Fitted Values", ylab = "Sqrt(|Residuals|)",
             main = "Scale-Location Plot")
        abline(h = 0, col = "red", lwd = 2)
        
        par(mfrow = c(1, 1))
        
        
        # Moran's I test on residuals 
        residuals <- residuals(gam_model)
        coords <- cbind(pf.env$longitude, pf.env$latitude)
        nb <- knn2nb(knearneigh(coords, k = 4))
        listw <- nb2listw(nb, style = "W")
        moran_test <- moran.test(residuals, listw)
          # Save Moran's I significance data 
        moran_p <- moran_test$p.value
        
        # Calculate the correlation matrix
        correlation_matrix <- cor(numeric_data, use = "complete.obs")
        
        # Identify all pairs with high correlation (absolute correlation > threshold)
        high_corr <- which(abs(correlation_matrix) > correlation_threshold & upper.tri(correlation_matrix), arr.ind = TRUE)
        
        # Create a data frame to store all correlated pairs and their correlation values
        correlated_pairs <- data.frame(
          predictor1 = rownames(correlation_matrix)[high_corr[, 1]],
          predictor2 = colnames(correlation_matrix)[high_corr[, 2]],
          correlation = correlation_matrix[high_corr]
        )
        
        # Print each high correlation pair for verification
        if(nrow(high_corr) !=0){
          for (i in 1:nrow(correlated_pairs)) {
            message <- sprintf("High correlation value of %.2f between %s and %s at edge %s and buf %s",
                               correlated_pairs$correlation[i],
                               correlated_pairs$predictor1[i],
                               correlated_pairs$predictor2[i], 
                               edge, 
                               buf)
            print(message)
          }
        }else{
          message <- sprintf("No highly correlated variables at edge %s and buf %s",
                             edge, 
                             buf)        
          # Now `correlated_pairs` contains all pairs with high correlation values for later use
          print(message)
        }
        
        # partial effects visualization ####
        # Clean up the interaction term to remove any spaces for accurate matching
        cleaned_interaction_term <- gsub(",\\s*", ",", cleaned_interaction_term)
        
        # Loop through each selected variable and assign p-values
        
        p_values <- sapply(selected_columns, function(var) {
          # Check if it's a single variable term in gam_summary$s.table
          if (paste0("s(", var, ")") %in% rownames(gam_summary$s.table)) {
            return(gam_summary$s.table[paste0("s(", var, ")"), 4])
          } 
          
          # If multiple interaction terms exist, loop through them and extract matching p-values
          else if (var == "man_lwr") {
           return( gam_summary$s.table[paste0("ti(edge_l.mangrove,land_water_ratio)"), 4])

            # If multiple p-values exist, return them as a concatenated string
          } else if  (var =="man_edge" ){
            return(gam_summary$s.table[paste0("ti(edge_l.mangrove,edge_perc)"), 4])
            
          } else {
            return(NA)  # Return NA if the variable or interaction term is not in the model
          }
        })
        

        # Loop through each selected variable to extract partial effects
        partial_selected_columns <- setdiff(selected_columns, c("man_lwr", "man_edge")) # not getting partial effects for interaction terms.
        for (var in partial_selected_columns) {
          print(paste("Processing variable:", var))
          index <- which(rownames(gam_summary$s.table) == paste0("s(", var, ")"))
          # if the variable is in the model and is significant
          if (length(index) > 0 && !is.na(p_values[var]) && p_values[var] < 0.05) {
            print(paste("Generating plot for:", var))
            
            # Get all variables except the one being varied
            other_vars <- setdiff(partial_selected_columns, var)
            
            # Create a sequence for the selected variable
            var_seq <- seq(min(pf.env[[var]], na.rm = TRUE), max(pf.env[[var]], na.rm = TRUE), length.out = 100)
            
            # If there are other variables, keep them constant at their mean
            if (length(other_vars) > 0) {
              pred_data  <- pf.env %>%
                summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE), .names = "{col}")) %>%
                slice(rep(1, 100)) %>%
                mutate(!!var := seq(min(pf.env[[var]], na.rm = TRUE),
                                    max(pf.env[[var]], na.rm = TRUE),
                                    length.out = 100)) %>%  # Vary only the selected predictor
                select(any_of(c(var, other_vars)))  # Keep only the relevant columns
            } else {
              pred_data <- data.frame(var_value = var_seq)
            }
            
            # Rename var_value column to match the variable name in the model
            colnames(pred_data)[colnames(pred_data) == "var_value"] <- var
            
            # Predict GAM response
            # Transform predictions back to original response scale
            pred_data$predicted_effect <- exp(predict(gam_model, newdata = pred_data, type = "link"))
            
            # Store for combined plotting
            pred_data$edge <- edge
            pred_data$buffer <- buf
            pred_data$variable <- var
            pred_data$r_sq <- r_sq
            pred_data$buffer_edge <- factor(paste0("Buffer: ", pred_data$buffer, " Edge:", pred_data$edge))
            effect_data <- rbind(effect_data, pred_data)
            
          }
        }
        ######
        
        
        
        # Dynamically create a new row for r.tab with selected columns and their values
        newline <- as.list(c(
          edge = edge,
          buffer = buf,
          sp = species_code,
          r_sq = r_sq,
          dev_ex = dev_ex,
          moran_p = moran_p,
          cov = max(correlated_pairs$correlation, na.rm = TRUE),
          setNames(as.list(p_values), paste0("p.", selected_columns))
        ))
        
        # Append the new row to the results table
        r.tab <- rbind(r.tab, newline, stringsAsFactors = FALSE)
      }, error = function(e) { cat("ERROR:", conditionMessage(e), "\n") })
    }
  }
  
  # Now that we've looped over every scale and stored all the data, we can plot and store it 
  
  # Define the output folder for the species
  output_folder <- file.path("./species_analysis", "output", "partial_effects", paste0(species_code, "_", scale))
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  # Define labels mapping
  labels_map <- c(
    "edge_perc" = '% Edge',
    "edge_l.mangrove" = "% Mangrove cover",
    "land_water_ratio" = "% Land"
  )
  
  # Partial effects plotting ########
  for(var in unique(effect_data$variable)){
    subset_effect <- subset(effect_data, variable == var)  # Check if var exists in the labels_map, otherwise use var itself
    x_label <- ifelse(var %in% names(labels_map), labels_map[[var]], var)
    x_label <- paste("z-scored", tolower(x_label))
    subset_effect$buffer_edge <- factor(subset_effect$buffer_edge, 
                                        levels = subset_effect %>%
                                          arrange(buffer, edge) %>%  # Sort by buffer first, then edge
                                          distinct(buffer_edge) %>%
                                          pull(buffer_edge))
    max_rsq = max(subset_effect$r_sq)
    # Plot all effects on the same graph with a gradient color scale
    if (nrow(subset_effect) > 0) {
      ggplot(subset_effect, aes(x = get(var), y = predicted_effect, group = buffer_edge, color = buffer)) +
        geom_line(size = 1.2) +
        # plot the 'best' Rsq model in bright green
        geom_line(data = subset(subset_effect, r_sq == max_rsq), 
                  aes(x = get(var), y = predicted_effect, linetype = "highest_r2"), 
                  color = "green", size = 1.2) + 
        scale_color_viridis_c(name = "Habitat radius (m)", option = "plasma") +  
        scale_linetype_manual(name = "", values = c("highest_r2" = "twodash"), 
                              labels = c("highest_r2" = "Highest R²\nscale")) +
        labs(x = x_label, y = "Partial effect") +
        theme_classic() +
        theme(text = element_text(size = 20), 
              legend.text = element_text(size = 17),  
              legend.title = element_text(size = 19, margin = margin(b = 10)),  # Add space below title
              legend.spacing.y = unit(0.5, "cm")  # Add vertical spacing in the legend
        )
      
      # Save the plot
      output_file <- file.path(output_folder, paste0("combined_partialeffect_", species_code, "_", var, ".png"))
      ggsave(output_file, width = 6, height = 6)

    }
  }
  
  # Save the data #########
  # Write the results to a CSV file
  write.csv(r.tab, file.path("./species_analysis", "output", "gam_stats", paste(scale,"_", species_code, "_outputs_over_scale.csv", sep = "")))
}
####### Application ######### 
species_list <- c("CALSAP", "PENSET", "PALSP", "CTEBOL") #define a list of species for the code to run through


# For loop for all the species in the list 
for(species_code in species_list){
  print(paste("Beginning", species_code, "analysis"))
  
  # Define the buffers and edges depending on which set scale you're running the analysis one. 
if(run_scale == "satscale"){
  buffers <- c(100, 150, 200, 250, 300, 400, 500, 600)# Buffer radii in meters
  edges <- c(1, 3, 5)# Edge distances in meters
  
}else{
  buffers <- c(20, 30, 50, 70, 100, 120, 150)  # Buffer radii in meters
  edges <- c(1, 3, 5)  # Edge distances in meters
}

# Load site data
pffw_sites <- read.csv("landscape_analysis/input/dropfield_2209.csv")

# Load community composition data and merge with environmental variables
comtab <- read.csv("./species_analysis/input/PtFouSept2022count_ns.csv") %>%
  replace(is.na(.), 0) %>%  # Replace NA values with 0
  filter(FUNGRA <= 30) %>%  # Filter out outlier site, sample taken when tide was too low
  merge(., pffw_sites[, c("site_date_key", "salinity", "temperature", "latitude", "longitude")], by = "site_date_key")

# Remove outliers dynamically
comtab <- comtab %>%
  filter(!(species_code %in% outlier_table$species & site_date_key %in% outlier_table$site_date_key))
# Run the model_compile function
results <- model_compile(comtab, 
                         species_code = species_code,
                         buffers = buffers,
                         edges = edges)

write.csv(results,  file.path("./species_analysis", "output", "gam_stats", "best_models", paste(run_scale, species_code,"model_list", ".csv", sep = "")))
          
# check results
print(paste(species_code, "MODEL_RESULTS"))
print(results)



# Count the number of times each model appears. 
unique_model_counts <- results %>%
  group_by(best_model) %>%
  summarise(count = n(), r.sq = mean(r.sq), mean_bufer = mean(buffer))

write.csv(unique_model_counts,  file.path("./species_analysis", "output", "gam_stats", "best_models", paste(run_scale, species_code,"model_comparison", ".csv", sep = "")))

# Display the counts of unique models
print(paste(species_code, "MODEL_COUNT"))
print(unique_model_counts)

# Select the model formula with the highest count (most frequent best model), if two models have the same count, the model with the highest Rsq is used. 
top_model <- unique_model_counts %>%
  arrange(desc(count),desc(r.sq)) %>%
  slice(1) %>%
  pull(best_model)

# Convert the model formula from text to an actual formula object
top_model_formula <- as.formula(top_model)
#Save top model in formula table
formula_table[[species_code]][[run_scale]] <- top_model_formula


print(paste("TOP FORMULA:",top_model))

write.table(paste(deparse(top_model_formula), collapse = " "), file.path("./species_analysis", "output", "gam_stats", paste(run_scale, species_code,"_topmodel", ".txt", sep = "")))

# Run the univariate_edge function with the top model
univariate_edge(comtab = comtab, 
                species_code = species_code,
                edges = edges, 
                buffers = buffers, 
                top_model_formula = top_model_formula)


# Read in the results from the GAM models
modelouts22 <- read.csv(file.path("./species_analysis", "output", "gam_stats", paste(run_scale,"_", species_code,"_outputs_over_scale", ".csv", sep = ""))) %>%
  mutate(species_code = species_code)
# List of variables to check for significance
variables <- c("cov", "p.edge_perc", "p.edge_l.mangrove", "p.man_lwr","p.man_edge", "p.land_water_ratio", "moran_p")
available_variables <- variables[variables %in% colnames(modelouts22)]

# Dynamically create the significance columns based on the presence of each variable
reg_modelouts <- modelouts22 %>%
  mutate(
    across(all_of(available_variables), 
           ~ ifelse(. <= 0.05, "y", "n"), 
           .names = "sig.{col}"),
    # Special condition for 'cov' column: check if cov is above 0.8
    sig.cov = ifelse(cov > 0.8, "y", "n")
  )

names(reg_modelouts) <- gsub("^sig\\.p\\.", "sig.", names(reg_modelouts))

# make sure breaks fit the scale
if(run_scale == "satscale"){
  breaks = seq(0, 600, by = 100)
} else {
  breaks = seq(0, 200, by = 10)
}
# Identify significant variables present in the dataset
sig_vars <- names(reg_modelouts)[grepl("^sig\\.", names(reg_modelouts)) & colSums(reg_modelouts == "y") > 0]

# Define colors and sizes dynamically
sig_colors <- setNames(rep("black", length(sig_vars)), sig_vars)
sig_sizes <- setNames(rep(4, length(sig_vars)), sig_vars)

# Manually adjust colors for sig.cov and sig.moran_p
if ("sig.cov" %in% sig_vars) sig_colors["sig.cov"] <- "red"
if ("sig.moran_p" %in% sig_vars) sig_colors["sig.moran_p"] <- "red"

# Adjust size for sig.cov
if ("sig.moran_p" %in% sig_vars) sig_sizes["sig.moran_p"] <- 2


# plot Rsq values over all scales. Use geom_point shapes to indicate significance of variables 
output_file <- file.path('./species_analysis/output/plots/', paste(run_scale, "_", species_code, '.png', sep = ""))

png(output_file, width = 14, height = 6, units = 'in', res = 300)

# Initialize ggplot object with no base data
p <- ggplot() +
  
  # Add lines for model results with R² by buffer, colored and linetyped by edge distance
  geom_line(
    reg_modelouts, 
    mapping = aes(x = buffer, y = r_sq, col = as.factor(edge), linetype = as.factor(edge)), 
    linewidth = 1.25
  ) +
  
  # Add points for significant land-water ratio
  geom_point(
    reg_modelouts[reg_modelouts$sig.land_water_ratio == "y",], 
    mapping = aes(x = buffer, y = r_sq, shape = "sig.land_water_ratio"), 
    col = "black", 
    size = 5
  ) +
  
  # Add points for significant interaction of % mangrove and land-water ratio
  geom_point(
    reg_modelouts[reg_modelouts$sig.man_lwr == "y",], 
    mapping = aes(x = buffer, y = r_sq, shape = "sig.man_lwr"), 
    col = "black", 
    size = 5
  ) +
  
  # Add points for significant % edge
  geom_point(
    reg_modelouts[reg_modelouts$sig.edge_perc == "y",], 
    mapping = aes(x = buffer, y = r_sq, shape = "sig.edge_perc"), 
    col = "black", 
    size = 5
  ) +
  
  # Add points for significant % mangrove cover
  geom_point(
    reg_modelouts[reg_modelouts$sig.edge_l.mangrove == "y",], 
    mapping = aes(x = buffer, y = r_sq, shape = "sig.edge_l.mangrove"), 
    size = 4
  ) +
  
  # Add points for significant interaction of % mangrove and % edge
  geom_point(
    reg_modelouts[reg_modelouts$sig.man_edge == "y",], 
    mapping = aes(x = buffer, y = r_sq, shape = "sig.man_edge"), 
    size = 4
  ) +
  
  # Add points for significant spatial covariance in residuals (Moran's I)
  geom_point(
    reg_modelouts[reg_modelouts$sig.moran_p == "y",], 
    mapping = aes(x = buffer, y = r_sq, shape = "sig.moran_p"), 
    size = 2, 
    col = 'red'
  ) +
  
  # Add points for significant covariance between % land and % edge
  geom_point(
    reg_modelouts[reg_modelouts$sig.cov == "y",], 
    mapping = aes(x = buffer, y = r_sq, shape = "sig.cov"), 
    col = "red", 
    size = 6
  ) +
  
  # Set axis labels with formatted text for R²
  labs(
    x = "Habitat radius (m)", 
    col = "Edge distance (m)", 
    y = expression("Adjusted GAM R"^2)
  ) + 
  
  # Customize x-axis with predefined breaks
  scale_x_continuous(breaks = breaks) +
  
  # Customize y-axis with fixed breaks and limits
  scale_y_continuous(breaks = seq(0, 0.7, by = 0.1), limits = c(-0.025,0.7)) +
  
  # Use grayscale color scale for lines (from dark to light)
  scale_color_grey(start = 0.2, end = 0.8) +
  
  # Manually define line types for edge distances
  scale_linetype_manual(
    values = c("1" = "dashed", "3" = "dotted", "5" = "dotdash")
  ) +
  
  # Customize legends: 
  #  - Adjust shape legend size/colors using external vectors `sig_sizes` & `sig_colors`
  #  - Override line legend linetypes for clarity
  guides(
    shape = guide_legend(title = "Significance", override.aes = list(size = as.numeric(sig_sizes), col = sig_colors)),
    col = guide_legend(override.aes = list(linetype = c("dashed", "dotted", "dotdash"))),
    linetype = "none"
  ) +
  
  # Manually assign shapes to each significance variable with clear labels
  scale_shape_manual(
    values = c("sig.edge_perc" = 5, 
               "sig.edge_l.mangrove" = 3, 
               "sig.moran_p" = 16,
               "sig.land_water_ratio" = 4,
               "sig.man_edge" = 2,
               "sig.man_lwr" = 2,
               "sig.cov" = 1),
    
    # Provide descriptive labels for each significance indicator
    labels = c(
      "sig.cov" = 'Covariance between % Land and % Edge',
      "sig.edge_perc" = '% Edge',
      "sig.edge_l.mangrove" = "% Mangrove cover",
      "sig.land_water_ratio" = "% Land",
      "sig.man_edge" = "Interaction of % mangrove and % edge",
      "sig.man_lwr" = "Interaction of % mangrove and % edge",
      "sig.moran_p" = "Spatial covariance in residuals"
    )
  ) +
  
  # Use a clean, minimal theme
  theme_classic() +
  
  # Adjust text size for readability
  theme(text = element_text(size = 18))


print(p)

dev.off()
}

#### SPECIES PLOTTING ######### 
# this will produce some species-specific plots depending on the species and scale you're plotting. 
# Plots are not consistent across the board, but differ based on needs for the manuscript 
for(species_code in species_list)
#### CTEBOL ####
if(species_code == "CTEBOL"){

if(run_scale == "satscale"){
  
  
  ##### land plot #####
  edge <- 1
  buf <- 400
  species_code <- "CTEBOL"
  run_scale = "satscale"
  # name scale habitat data file 
  habx <- paste("google2022_", "edge", edge, "_buf", buf, sep = "")
  
  # Load habitat data
  hab <-  read.csv(file.path("./landscape_analysis", "output", "satscale", paste(habx,".csv", sep = "")))
  hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
  
  # Merge habitat data with community composition table
  pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
  # Standardize variables and store statistics for back-transformation
  predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
  standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
  
  for (col in predictor_columns) {
    standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
    standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
    pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
  }
  # run the model selected above at the scale chosen and save the model
  best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
  summary(best_model)
  
  vis.gam(best_model)
  library(mgcv)
  
  # Calculate Cook's distance for the best_model
  cooks_distances <- cooks.distance(best_model)
  
  # Examine the Cook's distances
  summary(cooks_distances)
  
  # Plot Cook's distances
  plot(cooks_distances, type = "h", 
       main = "Cook's Distance", 
       ylab = "Cook's Distance", 
       xlab = "Observation Index")
  abline(h = 0.8, col = "red", lty = 2) # Threshold
  
  # Residual plot
  residuals_gam <- residuals(best_model, type = "response")
  fitted_values <- fitted(best_model)
  
  plot(fitted_values, residuals_gam, 
       xlab = "Fitted Values", ylab = "Residuals", 
       main = "Residuals vs Fitted Values")
  abline(h = 0, col = "red")
  
  
  #extract coefficients
  smooth_table<-summary(best_model)$s.table
  lwr22 = smooth_table["s(land_water_ratio)",4]
  edge22 <- smooth_table["s(edge_perc)",4]
  # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
  predict_tab_edg <- create_sequence_table(pf.env, "land_water_ratio") # variable in the function will be sequenced, all others in the datatable will be held at the mean
  
  # predict response values over the range of variable of interest
  predict_tab_edg$predicted_count <- predict(best_model, newdata = predict_tab_edg, type = "response")
  
  for (col in predictor_columns) {
    predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
    pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
  }
  
  
  png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'lwr.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
  
  # plot
  ggplot(predict_tab_edg, aes(x = land_water_ratio, y = predicted_count, col = edge_l.mangrove)) +
    geom_line(color = "darkgreen", size = 1) +
    geom_point(data = pf.env, aes(x = land_water_ratio, y = !!ensym(species_code)), size = 4) +
    annotate("text", x = mean(pf.env$land_water_ratio), y = max(pf.env[,as.character(species_code)]*0.8),
             label = paste("p =", round(lwr22, 4))) +
    labs(
      x = "% Land",
      y = "Darter goby abundance",
      col = "% Mangrove cover",
      title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                      edge, 
                      buf)        
    ) +
    scale_color_binned(type = 'viridis', labels = percent_format(scale = 100, suffix = "%", accuracy = 1))+
    scale_x_continuous(labels =percent_format(accuracy =.1)) +
    theme_classic() +
    theme(
      text = element_text(size = 16),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5)
    )
  
  
  dev.off()
  
  
}

if(run_scale == "smallscale"){
  ##### land plot #####
  edge <- 1 
  buf <- 150
  species_code <- "CTEBOL"
  run_scale = "smallscale"
  # name scale habitat data file 
  habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
  
  # Load habitat data
  hab <-  read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx,".csv", sep = "")))
  hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
  
  # Merge habitat data with community composition table
  pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
  # Standardize variables and store statistics for back-transformation
  predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
  standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
  
  for (col in predictor_columns) {
    standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
    standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
    pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
  }
  # run the model selected above at the scale chosen and save the model
  best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
  summary(best_model)
  # Calculate Cook's distance for the best_model
  cooks_distances <- cooks.distance(best_model)
  
  # Examine the Cook's distances
  summary(cooks_distances)
  
  # Plot Cook's distances
  plot(cooks_distances, type = "h", 
       main = "Cook's Distance", 
       ylab = "Cook's Distance", 
       xlab = "Observation Index")
  abline(h = 0.8, col = "red", lty = 2) # Threshold
  
  
  # Residual plot
  residuals_gam <- residuals(best_model, type = "response")
  fitted_values <- fitted(best_model)
  
  plot(fitted_values, residuals_gam, 
       xlab = "Fitted Values", ylab = "Residuals", 
       main = "Residuals vs Fitted Values")
  abline(h = 0, col = "red")
  
  
  #extract coefficients
  smooth_table<-summary(best_model)$s.table
  lwr22 = smooth_table["s(land_water_ratio)",4]
  # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
  predict_tab_edg <- create_sequence_table(pf.env, "land_water_ratio") # variable in the function will be sequenced, all others in the datatable will be held at the mean
  
  # predict response values over the range of variable of interest
  predict_tab_edg$predicted_count <- predict(best_model, newdata = predict_tab_edg, type = "response")
  
  for (col in predictor_columns) {
    predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
    pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
  }
  
  
  png(file.path('./species_analysis/output/plots/single_scale', paste('smallscale_', species_code,'lwr.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
  
  # plot
  ggplot(predict_tab_edg, aes(x = land_water_ratio, y = predicted_count, col = edge_l.mangrove)) +
    geom_line(color = "darkgreen", size = 1) +
    geom_point(data = pf.env, aes(x = land_water_ratio, y = !!ensym(species_code)), size = 4) +
    annotate("text", x = mean(pf.env$land_water_ratio), y = max(pf.env[,as.character(species_code)]*0.8),
             label = paste("p =", round(lwr22, 3))) +
    labs(
      x = "% Land",
      y = "Darter goby abundance",
      col = "% Mangrove cover",
      title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                      edge, 
                      buf)        
    ) +
    scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
    scale_x_continuous(labels = percent_format(accuracy = 1))+
    theme_classic() +
    theme(
      text = element_text(size = 16),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5)
    )
  
  
  dev.off()
  
  
} 
  
}
#### PENSET ####
if(species_code == "PENSET"){
  
  if(run_scale == "smallscale"){
    ##### mangrove plot ###############
    
    # Based on GAM outputs, select a scale to plot species abundance
    edge <- 1 
    buf <- 100
    species_code <- "PENSET"
    run_scale = "smallscale"
    # name scale habitat data file 
    habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    
    vis.gam(best_model, theta = 120, n.grid = 50, lwd = 0.4)
    
    
    # Calculate Cook's distance for the best_model
    cooks_distances <- cooks.distance(best_model)
    
    # Examine the Cook's distances
    summary(cooks_distances)
    
    # Plot Cook's distances
    plot(cooks_distances, type = "h", 
         main = "Cook's Distance", 
         ylab = "Cook's Distance", 
         xlab = "Observation Index")
    abline(h = 4 / length(cooks_distances), col = "red", lty = 2) # Threshold
    
    threshold <- 4 / nrow(pf.env)
    
    # Identify influential observations
    influential_obs <- which(cooks_distances > threshold)
    
    # Print influential observations
    print(influential_obs)
    weights <- rep(1, nrow(pf.env)) # ADDED DUE TO HIGH OUTLIERS at obs 1 and 49 (cooks distance of 0.4 and 0.86 respectively)
   # weights[c(1)] <- 0.1  # Assign lower weight to influential observation
    
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    man22<-smooth_table["s(edge_l.mangrove)",4]
    
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab <- create_sequence_table(pf.env, "edge_l.mangrove") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    # predict response values over the range of variable of interest
    predict_tab$predicted_count <- predict(best_model, newdata = predict_tab, type = "response")
    
    # Back-transform each variable in predict_tab
    for (col in predictor_columns) {
      predict_tab[[col]] <- back_transform(predict_tab[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    # Generate predictions on the original data (pf.env)
    pf.env$predicted_count <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    library(ggplot2)
    
    # Plot observed vs. predicted values
    ggplot(data = pf.env, aes(x = .data[[species_code]], y = predicted_count)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, color = "red") +
      labs(x = "Observed", y = "Predicted") +
      theme_classic()
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('smallscale_', species_code,'mangrove.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab, aes(x = edge_l.mangrove, y = predicted_count)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = edge_l.mangrove, y = !!ensym(species_code), col = (edge_perc)*1000), size = 4) +
      annotate("text", x = mean(pf.env$edge_l.mangrove) * 0.5, y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(man22, 3))) +
      labs(
        x = "% Mangrove edge",
        y = "White shrimp abundance",
        color = "‰ Edge",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = label_number(suffix = "‰")) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    ##### edge plot #####
    edge <- 1 
    buf <- 30
    species_code <- "PENSET"
    
    # name scale habitat data file 
    habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    
    # Calculate Cook's distance for the best_model
    cooks_distances <- cooks.distance(best_model)
    
    # Examine the Cook's distances
    summary(cooks_distances)
    
    # Plot Cook's distances
    plot(cooks_distances, type = "h", 
         main = "Cook's Distance", 
         ylab = "Cook's Distance", 
         xlab = "Observation Index")
    abline(h = 0.6, col = "red", lty = 2) # Threshold
    
    threshold <- 4 / nrow(pf.env)
    
    # Identify influential observations
    influential_obs <- which(cooks_distances > threshold)
    
    # Print influential observations
    print(influential_obs)
    weights <- rep(1, nrow(pf.env)) # ADDED DUE TO HIGH OUTLIERS at obs 1 and 49 (cooks distance of 0.4 and 0.86 respectively)
 #   weights[c(1)] <- 0.1  # Assign lower weight to influential observation
    
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    man22<-smooth_table["s(edge_l.mangrove)",4]
    edge22 = smooth_table["s(edge_perc)", 4]
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab_edg <- create_sequence_table(pf.env, "edge_perc") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    
    # predict response values over the range of variable of interest
    predict_tab_edg$predicted_count <- predict(best_model, newdata = predict_tab_edg, type = "response")

    for (col in predictor_columns) {
      predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    pf.env$full_predict<- predict(best_model, newdata = pf.env, type = "response")%>% rollmean( k = 5, fill = NA)
    png(file.path('./species_analysis/output/plots/single_scale', paste('smallscale_', species_code,'edge.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab_edg, aes(x = edge_perc*1000, y = predicted_count, col = edge_l.mangrove)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_line(data = pf.env, aes(y = full_predict, x = edge_perc*1000), color = ) +  
      geom_point(data = pf.env, aes(x = edge_perc*1000, y = !!ensym(species_code)), size = 4) +
      annotate("text", x = mean(pf.env$edge_perc*1000) * 0.5, y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(edge22, 4))) +
      labs(
        x = "‰ Edge",
        y = "White shrimp abundance",
        col = "% Mangrove cover",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      scale_x_continuous(labels = label_number(suffix = "‰")) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    
    
  }
  if(run_scale == "satscale"){
    #### mangrove ####
    # Based on GAM outputs, select a scale to plot species abundance
    edge <- 1
    buf <- 300
    species_code <- "PENSET"
    
    # name scale habitat data file 
    habx <- paste("google2022_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "satscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    
    vis.gam(best_model)
    # Calculate Cook's distance for the best_model
    cooks_distances <- cooks.distance(best_model)
    
    # Examine the Cook's distances
    summary(cooks_distances)
    
    # Plot Cook's distances
    plot(cooks_distances, type = "h", 
         main = "Cook's Distance", 
         ylab = "Cook's Distance", 
         xlab = "Observation Index")
    abline(h = 0.8, col = "red", lty = 2) # Threshold
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    man22<-smooth_table["s(edge_l.mangrove)",4]
    
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab <- create_sequence_table(pf.env, "edge_l.mangrove") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    # predict response values over the range of variable of interest
    predict_tab$predicted_count <- predict(best_model, newdata = predict_tab, type = "response")
    
    # Back-transform each variable in predict_tab
    for (col in predictor_columns) {
      predict_tab[[col]] <- back_transform(predict_tab[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    # Generate predictions on the original data (pf.env)
    pf.env$predicted_count <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    library(ggplot2)
    
    # Plot observed vs. predicted values
    ggplot(data = pf.env, aes(x = .data[[species_code]], y = predicted_count)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, color = "red") +
      labs(x = "Observed", y = "Predicted") +
      theme_classic()
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'mangrove.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab, aes(x = edge_l.mangrove, y = predicted_count, col = edge_perc*1000)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = edge_l.mangrove, y = !!ensym(species_code), col = edge_perc*1000), size = 4) +
      annotate("text", x = mean(pf.env$edge_l.mangrove) * 0.5, y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(man22, 3))) +
      labs(
        x = "% Mangrove edge",
        y = "White shrimp abundance",
        color = "‰ Edge",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = label_number(suffix = "‰")) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    ##### edge plot #####
    edge <- 1
    buf <- 300
    species_code <- "PENSET"
    
    # name scale habitat data file 
    habx <- paste("google2022_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "satscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    # Calculate Cook's distance for the best_model
    cooks_distances <- cooks.distance(best_model)
    
    # Examine the Cook's distances
    summary(cooks_distances)
    
    # Plot Cook's distances
    plot(cooks_distances, type = "h", 
         main = "Cook's Distance", 
         ylab = "Cook's Distance", 
         xlab = "Observation Index")
    abline(h = 4 / length(cooks_distances), col = "red", lty = 2) # Threshold
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    edge22 <- smooth_table["s(edge_perc)",4]
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab_edg <- create_sequence_table(pf.env, "edge_perc") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    
    # predict response values over the range of variable of interest
    predict_tab_edg$predicted_count <- predict(best_model, newdata = predict_tab_edg, type = "response")
    
 
    
    pf.env$true_predict <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    # Compute residuals
    pf.env$residuals <- pf.env$PENSET - pf.env$true_predict
    
    # Apply LOWESS smoothing
    loess_fit <- loess(true_predict ~ edge_perc, data = pf.env, span = 0.5)    
    pf.env$Smoothed_Pred <- predict(loess_fit)
    # Plot observed points colored by residuals and smoothed prediction curve
       for (col in predictor_columns) {
      predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'prediction.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
   
    ggplot(pf.env, aes(x = edge_perc, y = PENSET)) +
      geom_point(aes(color = edge_l.mangrove), size = 6, alpha = 1) + 
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      
      # Map linetype inside aes() to include Smoothed_Pred in the legend
      geom_line(aes(y = Smoothed_Pred, linetype = "Smoothed \nprediction"), color = "red", size = 1.2) +
      
      # Use scale_linetype_manual to ensure it appears in the legend
      scale_linetype_manual(name = "", values = c("Smoothed \nprediction" = "solid")) +
      
      labs(x = "‰ Edge", y = "White shrimp abundance", col = "% Mangrove cover",
           title = sprintf("Edge distance = %s m \n Buffer radius = %s m", edge, buf)) +
      scale_x_continuous(labels = percent_format(scale = 1000, suffix = "‰")) +
      
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    dev.off()
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'edge.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab_edg, aes(x = edge_perc, y = predicted_count, col = edge_l.mangrove)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = edge_perc, y = !!ensym(species_code)), size = 4) +
      annotate("text", x = mean(pf.env$edge_perc), y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(edge22, 6))) +
      labs(
        x = "‰ Edge",
        y = "White Shrimp abundance",
        col = "% Mangrove edge",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      scale_x_continuous(labels = percent_format(scale = 1000 ,suffix = "‰")) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
  }
}
#### CALSAP ####
if(species_code == "CALSAP"){
  if(run_scale == "satscale"){
    
    ##### edge plot #####
    edge <- 5
    buf <- 200
    species_code <- "CALSAP"
    
    # name scale habitat data file 
    habx <- paste("google2022_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "satscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    
    vis.gam(best_model)
    library(mgcv)
    
    # Calculate Cook's distance for the best_model
    cooks_distances <- cooks.distance(best_model)
    
    # Examine the Cook's distances
    #summary(cooks_distances)
    
    # Plot Cook's distances
    plot(cooks_distances, type = "h", 
      main = "Cook's Distance", 
       ylab = "Cook's Distance", 
       xlab = "Observation Index",
    h = 0.8, col = "red", lty = 2) # Threshold
    
    #threshold <- 0.8
    
    # Identify influential observations
    #influential_obs <- which(cooks_distances > threshold)
    
    # Print influential observations
    #print(influential_obs)
    #weights <- rep(1, nrow(pf.env)) # ADDED DUE TO HIGH OUTLIERS at obs 1 and 49 (cooks distance of 0.4 and 0.86 respectively)
    #[c(1)] <- 0.1  # Assign lower weight to influential observation
    
    #best_model <- gam(as.formula(top_model_formula), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML", weights = weights)
    #summary(best_model)
    
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    lwr22 = smooth_table["s(land_water_ratio)",4]
    edge22 <- smooth_table["s(edge_perc)",4]
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab_edg <- create_sequence_table(pf.env, "edge_perc") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    
    # predict response values over the range of variable of interest
    predict_tab_edg$predicted_count <- predict(best_model, newdata = predict_tab_edg, type = "response")
    
    for (col in predictor_columns) {
      predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    
    pf.env$true_predict <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    # Compute residuals
    pf.env$residuals <- pf.env$CALSAP - pf.env$true_predict
    
    # Apply LOWESS smoothing
    loess_fit <- loess(true_predict ~ land_water_ratio, data = pf.env, span = 0.5)    
    pf.env$Smoothed_Pred <- predict(loess_fit)
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'prediction.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    ggplot(pf.env, aes(x = land_water_ratio, y = CALSAP)) +
      geom_point(aes(color = edge_l.mangrove), size = 6, alpha = 1) + 
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      geom_line(aes(y = Smoothed_Pred), color = "red", size = 1.2) +
      labs(x = "‰ Edge", y = "Blue crab abundance", col = "% Mangrove cover") +
      scale_x_continuous(labels = percent_format(scale = 1000 ,suffix = "‰")) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    dev.off()
    
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'edge.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab_edg, aes(x = edge_perc*1000, y = predicted_count, col = land_water_ratio)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = edge_perc*1000, y = !!ensym(species_code)), size = 4) +
      annotate("text", x = mean(pf.env$edge_perc*1000), y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(edge22, 3))) +
      labs(
        x = "‰ Edge",
        y = "Blue crab abundance",
        col = "% Land",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      scale_x_continuous(labels = label_number(suffix = "‰")) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    ##### land plot #####
    edge <- 5
    buf <- 150
    species_code <- "CALSAP"
    run_scale == "satscale"
    # name scale habitat data file 
    habx <- paste("google2022_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "satscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    
    vis.gam(best_model)
    library(mgcv)
    
    # Calculate Cook's distance for the best_model
    cooks_distances <- cooks.distance(best_model)
    
    # Examine the Cook's distances
    summary(cooks_distances)
    
    # Plot Cook's distances
    plot(cooks_distances, type = "h", 
         main = "Cook's Distance", 
         ylab = "Cook's Distance", 
         xlab = "Observation Index")
    abline(h = 0.5 , col = "red", lty = 2) # Threshold
    
    #threshold <- 0.8
    
    # Identify influential observations
    #influential_obs <- which(cooks_distances > threshold)
    
    # Print influential observations
    #print(influential_obs)
    #weights <- rep(1, nrow(pf.env)) # ADDED DUE TO HIGH OUTLIERS at obs 1 and 49 (cooks distance of 0.4 and 0.86 respectively)
    #[c(1)] <- 0.1  # Assign lower weight to influential observation
    
    #best_model <- gam(as.formula(top_model_formula), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML", weights = weights)
    #summary(best_model)
    
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    lwr22 = smooth_table["s(land_water_ratio)",4]
    edge22 <- smooth_table["s(edge_perc)",4]
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab_edg <- create_sequence_table(pf.env, "land_water_ratio") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    
    # predict response values over the range of variable of interest
    predict_tab_edg$predicted_count <- predict(best_model, newdata = predict_tab_edg, type = "response")
   
    
    pf.env$true_predict <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    # Compute residuals
    pf.env$residuals <- pf.env$CALSAP - pf.env$true_predict
    
    # Apply LOWESS smoothing
    loess_fit <- loess(true_predict ~ land_water_ratio, data = pf.env, span = 0.5)    
    pf.env$Smoothed_Pred <- predict(loess_fit)
    # Plot observed points colored by residuals and smoothed prediction curve
    
    for (col in predictor_columns) {
      predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'prediction.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    ggplot(pf.env, aes(x = land_water_ratio, y = CALSAP)) +
      geom_point(aes(color = edge_perc), size = 6, alpha = 1) + 
      scale_color_binned(type = 'viridis', labels = percent_format(scale = 1000, suffix = "‰")) +
      
      # Map linetype inside aes() to include Smoothed_Pred in the legend
      geom_line(aes(y = Smoothed_Pred, linetype = "Smoothed \nprediction"), color = "red", size = 1.2) +
      
      # Use scale_linetype_manual to ensure it appears in the legend
      scale_linetype_manual(name = "", values = c("Smoothed \nprediction" = "solid")) +
      
      labs(x = "% Land", y = "Blue crab abundance", col = "‰ Edge",
           title = sprintf("Edge distance = %s m \n Buffer radius = %s m", edge, buf)) +
      scale_x_continuous(labels = percent_format()) +

      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    dev.off()
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'lwr.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab_edg, aes(x = land_water_ratio, y = predicted_count, col = edge_perc)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = land_water_ratio, y = !!ensym(species_code)), size = 4) +
      annotate("text", x = mean(pf.env$land_water_ratio), y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(lwr22, 3))) +
      labs(
        x = "% Land",
        y = "Blue crab abundance",
        col = "‰ Edge",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(scale = 1000, suffix = "‰", accuracy = 0.001))+
      scale_x_continuous(labels =percent_format(accuracy =.1)) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
  }
  if(run_scale == "smallscale"){
    
    
    ########## mangrove ###############
    
    # Based on GAM outputs, select a scale to plot species abundance
    edge <- 3 
    buf <- 50
    species_code <- "CALSAP"
    
    # name scale habitat data file 
    habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(top_model_formula), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    
    
    
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    man22<-smooth_table["s(edge_l.mangrove)",4]
    
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab <- create_sequence_table(pf.env, "edge_l.mangrove") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    # predict response values over the range of variable of interest
    predict_tab$predicted_count <- predict(best_model, newdata = predict_tab, type = "response")
    
    # Back-transform each variable in predict_tab
    for (col in predictor_columns) {
      predict_tab[[col]] <- back_transform(predict_tab[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    # Generate predictions on the original data (pf.env)
    pf.env$predicted_count <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    library(ggplot2)
    
    # Plot observed vs. predicted values
    ggplot(data = pf.env, aes(x = .data[[species_code]], y = predicted_count)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, color = "red") +
      labs(x = "Observed", y = "Predicted") +
      theme_classic()
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('smallscale_', species_code,'mangrove.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab, aes(x = edge_l.mangrove, y = predicted_count)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = edge_l.mangrove, y = !!ensym(species_code), col = land_water_ratio), size = 4) +
      annotate("text", x = mean(pf.env$edge_l.mangrove) * 0.5, y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(man22, 3))) +
      labs(
        x = "% Mangrove cover",
        y = "Blue crab abundance",
        color = "‰ Land",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 0.5)) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    ##### land plot #####
    edge <- 5
    buf <- 150
    species_code <- "CALSAP"
    
    # name scale habitat data file 
    habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    # Calculate Cook's distance for the best_model
    cooks_distances <- cooks.distance(best_model)
    
    # Examine the Cook's distances
    summary(cooks_distances)
    
    # Plot Cook's distances
    plot(cooks_distances, type = "h", 
         main = "Cook's Distance", 
         ylab = "Cook's Distance", 
         xlab = "Observation Index")
    abline(h = 0.8, col = "red", lty = 2) # Threshold
    
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    lwr22 = smooth_table["s(land_water_ratio)",4]
    edge22 <- smooth_table["s(edge_perc)",4]
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab_edg <- create_sequence_table(pf.env, "land_water_ratio") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    
    # predict response values over the range of variable of interest
    predict_tab_edg$predicted_count <- predict(best_model, newdata = predict_tab_edg, type = "response")

    
    
    pf.env$true_predict <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    # Compute residuals
    pf.env$residuals <- pf.env$CALSAP - pf.env$true_predict
    
    # Apply LOWESS smoothing
    loess_fit <- loess(true_predict ~ land_water_ratio, data = pf.env, span = 0.5)    
    pf.env$Smoothed_Pred <- predict(loess_fit)
       
    for (col in predictor_columns) {
      predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
     
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'prediction.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    ggplot(pf.env, aes(x = land_water_ratio, y = CALSAP)) +
      geom_point(aes(color = edge_l.mangrove), size = 6, alpha = 1) + 
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      geom_line(aes(y = Smoothed_Pred), color = "red", size = 1.2) +
      labs(x = "% Edge", y = "Blue crab abundance", col = "% Mangrove cover") +
      scale_x_continuous(labels = percent_format()) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    dev.off()
    png(file.path('./species_analysis/output/plots/single_scale', paste('smallscale_', species_code,'lwr_edge5_buf150.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab_edg, aes(x = land_water_ratio, y = predicted_count, col = edge_l.mangrove)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = land_water_ratio, y = !!ensym(species_code)), size = 4) +
      annotate("text", x = mean(pf.env$land_water_ratio), y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(lwr22, 4))) +
      labs(
        x = "% Land",
        y = "Blue crab abundance",
        col = "% Mangrove cover",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      scale_x_continuous(labels = percent_format(accuracy = 1))+
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    ##### land plot 2 #####
    edge <- 3
    buf <- 50
    species_code <- "CALSAP"
    run_scale = "smallscale"
    # name scale habitat data file 
    habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    # Calculate Cook's distance for the best_model
    cooks_distances <- cooks.distance(best_model)
    
    # Examine the Cook's distances
    summary(cooks_distances)
    
    # Plot Cook's distances
    plot(cooks_distances, type = "h", 
         main = "Cook's Distance", 
         ylab = "Cook's Distance", 
         xlab = "Observation Index")
    abline(h = 0.8, col = "red", lty = 2) # Threshold
    
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    lwr22 = smooth_table["s(land_water_ratio)",4]
    edge22 <- smooth_table["s(edge_perc)",4]
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab_edg <- create_sequence_table(pf.env, "land_water_ratio") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    
    # predict response values over the range of variable of interest
    predict_tab_edg$predicted_count <- predict(best_model, newdata = predict_tab_edg, type = "response")
    
    for (col in predictor_columns) {
      predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    
    
    pf.env$true_predict <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    # Compute residuals
    pf.env$residuals <- pf.env$CALSAP - pf.env$true_predict
    
    # Apply LOWESS smoothing
    loess_fit <- loess(true_predict ~ land_water_ratio, data = pf.env, span = 0.5)    
    pf.env$Smoothed_Pred <- predict(loess_fit)
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'prediction.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    ggplot(pf.env, aes(x = land_water_ratio, y = CALSAP)) +
      geom_point(aes(color = edge_l.mangrove), size = 6, alpha = 1) + 
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      geom_line(aes(y = Smoothed_Pred), color = "red", size = 1.2) +
      labs(x = "% Edge", y = "Blue crab abundance", col = "% Mangrove cover") +
      scale_x_continuous(labels = percent_format(scale = 1000 ,suffix = "‰")) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    dev.off()
    png(file.path('./species_analysis/output/plots/single_scale', paste('smallscale_', species_code,'lwr_edge3_buf50.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab_edg, aes(x = land_water_ratio, y = predicted_count, col = edge_l.mangrove)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = land_water_ratio, y = !!ensym(species_code)), size = 4) +
      annotate("text", x = mean(pf.env$land_water_ratio), y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(lwr22, 3))) +
      labs(
        x = "% Land",
        y = "Blue crab Abundance",
        col = "% Mangrove cover",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      scale_x_continuous(labels = percent_format(accuracy = 1))+
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    
    
    ##### edge plot #####
    edge <- 3
    buf <- 50
    species_code <- "CALSAP"
    
    # name scale habitat data file 
    habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(top_model_formula), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    edge22 <- smooth_table["s(edge_perc)",4]
    lwr22 = smooth_table["s(land_water_ratio)",4]
    
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab_edg <- create_sequence_table(pf.env, "edge_perc") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    
    # predict response values over the range of variable of interest
    predict_tab_edg$predicted_count <- predict(best_model, newdata = predict_tab_edg, type = "response")
    
    for (col in predictor_columns) {
      predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('smallscale_', species_code,'edge_edge3_buf50.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab_edg, aes(x = edge_perc*1000, y = predicted_count, col = land_water_ratio)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = edge_perc*1000, y = !!ensym(species_code)), size = 4) +
      annotate("text", x = mean(pf.env$edge_perc*1000), y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(edge22, 5))) +
      labs(
        x = "‰ Edge",
        y = "Blue crab abundance",
        col = "% Land",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      scale_x_continuous(labels = label_number(suffix = "‰")) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    
    
    ######## MULTISCALE ANALYSIS #########
    # r50 
    
    buf<-50
    edge <-3 
    
    habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
    hab <-  read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    r50_hab <- hab %>%
      rename_with(~paste0("r50_", .), everything())
    
    # r150 
    
    buf<-150
    edge <-3 
    
    habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
    hab <-  read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    r150_hab <- hab %>%
      rename_with(~paste0("r150_", .), everything())
    
    multiscale <- merge(r150_hab, r50_hab, by.x = "r150_site_date_key", by.y = "r50_site_date_key") %>% 
      merge( comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by.x = 'r150_site_date_key', by.y = "site_date_key") 
    
    all_predictors <- c("r50_edge_perc", "r50_edge_l.mangrove", 
                        "r50_land_water_ratio", "r50_edge_l.mangrove", 
                        "r150_edge_perc", "r150_edge_l.mangrove","r150_land_water_ratio")
    # Iterate over all subsets of predictors
    
    best_aic <- Inf
    models_within_aic <- list()  # Store models within the AIC threshold
    all_models_info <- data.frame(formula = character(), aic = numeric(), stringsAsFactors = FALSE)  # Store all tested models
    
    
    for (subset_size in 1:length(all_predictors)) {
      predictor_subsets <- combn(all_predictors, 2, simplify = FALSE)
      for (subset_predictors in predictor_subsets) {
        # Skip subsets containing highly correlated pairs
        # Create formula terms with smooth terms for each predictor
        formula_terms <- paste(sprintf("s(%s, k = 4)", subset_predictors), collapse = " + ")
        
        # Add interaction terms if applicable
        if ("r50_edge_l.mangrove" %in% subset_predictors) {
          if ("r150_land_water_ratio" %in% subset_predictors) {
            # Interaction between edge length of mangrove and land-water ratio
            formula_terms <- paste(formula_terms, " + ti(r50_edge_l.mangrove, r150_land_water_ratio, k = 4)")
          } else if ("r150_edge_perc" %in% subset_predictors) {
            # Interaction between edge length of mangrove and edge percentage
            formula_terms <- paste(formula_terms, " + ti(r50_edge_l.mangrove, r150_edge_perc, k = 4)")
          }
        }else if ("r150_edge_l.mangrove" %in% subset_predictors) {
          if ("r50_land_water_ratio" %in% subset_predictors) {
            # Interaction between edge length of mangrove and land-water ratio
            formula_terms <- paste(formula_terms, " + ti(r150_edge_l.mangrove, r50_land_water_ratio, k = 4)")
          } else if ("r50_edge_perc" %in% subset_predictors) {
            # Interaction between edge length of mangrove and edge percentage
            formula_terms <- paste(formula_terms, " + ti(r150_edge_l.mangrove, r50_edge_perc, k = 4)")
          }
        }
        
        # Construct the full formula and fit the GAM model
        formula <- as.formula(paste("CALSAP", "~", formula_terms))
        model <- gam(formula, data = multiscale, family = tw(link = "log"), method = "ML")
        model_aic <- AIC(model)
        aic_threshold <- 2
        
        # Store the model formula and AIC
        all_models_info <- rbind(all_models_info, data.frame(formula = paste(deparse(formula), collapse = " "), aic = model_aic, stringsAsFactors = FALSE))
        
        # Update the list of best models based on AIC
        if (model_aic < best_aic) {
          best_aic <- model_aic
          models_within_aic <- list(model)  # Reset list with the new best model
        } else if (model_aic <= best_aic + aic_threshold) {
          models_within_aic <- c(models_within_aic, list(model))  # Add model within AIC range
        }
      }
    }
    r.tab <- data.frame(
      edge = numeric(),
      buffer = numeric(),
      sp = character(),
      best_model = character(),
      aic = numeric(),
      stringsAsFactors = FALSE
    )
    for (i in seq_along(models_within_aic)) {
      model_info <- list(
        formula = formula(models_within_aic[[i]]),
        aic = AIC(models_within_aic[[i]])
      )
      
      # Create a new row with buffer, edge, species, model formula, and AIC
      newline <- data.frame(
        edge = edge,
        buffer = buf,
        sp = species_code,
        best_model = paste(deparse(model_info$formula), collapse = " "),
        aic = model_info$aic,
        stringsAsFactors = FALSE
      )
      
      # Append the new row to the results table
      r.tab <- rbind(r.tab, newline)
    }
    
    
    
    # run the model selected above at the scale chosen and save the model
    full_model <- gam(CALSAP ~ s(r50_edge_perc, k = 4) + s(r50_edge_l.mangrove, k = 4) + s(r50_land_water_ratio,k = 4) + ti(r50_edge_l.mangrove, r50_land_water_ratio, k = 4)+
                        s(r150_edge_perc, k = 4) + s(r150_edge_l.mangrove, k = 4) + s(r150_land_water_ratio,k = 4) + ti(r150_edge_l.mangrove, r150_land_water_ratio, k = 4)
                      , data = multiscale, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(full_model)
    
    best_model1 <-  gam(CALSAP ~ s(r50_edge_perc, k = 4)+ s(r150_land_water_ratio,k = 4) + ti(r150_edge_l.mangrove, r150_land_water_ratio, k = 4), 
                        data = multiscale, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model1)
    AIC(best_model1)
    AIC(full_model)
    
    plot(best_model1)
    ##### MUTLISCALE EDGEPERC ###########
    
    
    
    #where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab_edg <- create_sequence_table(multiscale, "r150_land_water_ratio") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    
    # predict response values over the range of variable of interest
    multiscale$predicted_count <- predict(best_model1, newdata = multiscale, type = "response")
    
    write.csv(multiscale,"CALSAP_multiscale.csv")
    
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('smallscale_', species_code,'MULTISCALE.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    
    
    # plot
    ggplot(multiscale, aes(x = r150_land_water_ratio, y = predicted_count, col = r50_edge_perc*1000)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = multiscale, aes(x = r150_land_water_ratio, y = !!ensym(species_code)), size = 4) +
      annotate("text", x = mean(multiscale$r150_land_water_ratio), y = max(multiscale[,as.character(species_code)]*0.8),
               label = paste("p =", round(edge22, 5))) +
      labs(
        col = "‰ Edge at 50m ",
        y = "Blue crab abundance",
        x = "% Land at 150m",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      scale_x_continuous(labels = label_number(suffix = "‰")) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    
    # plot
    ggplot(predict_tab_edg, aes(x = r50_edge_perc*1000, y = predicted_count, col = r150_land_water_ratio)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = multiscale, aes(x = r50_edge_perc*1000, y = !!ensym(species_code)), size = 4) +
      annotate("text", x = mean(multiscale$r50_edge_perc*1000), y = max(multiscale[,as.character(species_code)]*0.8),
               label = paste("p =", round(edge22, 5))) +
      labs(
        x = "‰ Edge at 50m ",
        y = "Blue crab abundance",
        col = "% Land at 150m",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      scale_x_continuous(labels = label_number(suffix = "‰")) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
  }
}
  
#### PALSP ####
if(species_code == "PALSP"){
  if(run_scale == "satscale"){
    
    #### mangrove ####
    # Based on GAM outputs, select a scale to plot species abundance
    edge <- 3
    buf <- 300
    species_code <- "PALSP"
    run_scale = "satscale"
    # name scale habitat data file 
    habx <- paste("google2022_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "satscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    
    vis.gam(best_model)
    library(mgcv)
    
    # Calculate Cook's distance for the best_model
    cooks_distances <- cooks.distance(best_model)
    
    # Examine the Cook's distances
    summary(cooks_distances)
    
    # Plot Cook's distances
    plot(cooks_distances, type = "h", 
         main = "Cook's Distance", 
         ylab = "Cook's Distance", 
         xlab = "Observation Index")
    abline(h = 4 / length(cooks_distances), col = "red", lty = 2) # Threshold
    
    threshold <- 4 / nrow(pf.env)
    
    # Identify influential observations
    influential_obs <- which(cooks_distances > threshold)
    
    # Print influential observations
    print(influential_obs)
    weights <- rep(1, nrow(pf.env)) # ADDED DUE TO HIGH OUTLIERS at obs 1 and 49 (cooks distance of 0.4 and 0.86 respectively)
    weights[c(1)] <- 0.01 # Assign lower weight to influential observation
    
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML", weights = weights)
    summary(best_model)
    # INITIAL RUN OF THIS ANALYSIS INDICATED A COOKS DISTANCE OF 0.86 for site_date_key 05244dc3, this value was down_weighhted in the final analysis for thios species
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    man22<-smooth_table["s(edge_l.mangrove)",4]
    
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab <- create_sequence_table(pf.env, "edge_l.mangrove") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    # predict response values over the range of variable of interest
    predict_tab$predicted_count <- predict(best_model, newdata = predict_tab, type = "response")
    
    # Back-transform each variable in predict_tab

    # Generate predictions on the original data (pf.env)
    pf.env$predicted_count <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    library(ggplot2)
    
    
    pf.env$true_predict <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    # Compute residuals
    pf.env$residuals <- pf.env$PALSP - pf.env$true_predict
    
    # Apply LOWESS smoothing
    loess_fit <- loess(true_predict ~ edge_l.mangrove, data = pf.env, span = 0.5)    
    pf.env$Smoothed_Pred <- predict(loess_fit)
    
    for (col in predictor_columns) {
      predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'prediction.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    ggplot(pf.env, aes(x = edge_l.mangrove, y = PALSP)) +
      geom_point(aes(color = edge_perc), size = 6, alpha = 1) + 
      scale_color_binned(type = 'viridis', labels = percent_format(scale = 1000, suffix = "‰")) +
      geom_line(aes(y = Smoothed_Pred), color = "red", size = 1.2) +
      labs(x = "% Edge", y = "White shrimp abundance", col = "% Mangrove cover") +
      scale_x_continuous(labels = percent_format()) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    dev.off()
    png(file.path('./species_analysis/output/plots/single_scale', paste('smallscale_', species_code,'lwr_edge5_buf150.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # Plot observed vs. predicted values
    ggplot(data = pf.env, aes(x = .data[[species_code]], y = predicted_count)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, color = "red") +
      labs(x = "Observed", y = "Predicted") +
      theme_classic()
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'mangrove.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab, aes(x = edge_l.mangrove, y = predicted_count)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = edge_l.mangrove, y = !!ensym(species_code), col = land_water_ratio), size = 4) +
      annotate("text", x = mean(pf.env$edge_l.mangrove) * 0.5, y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(man22, 4))) +
      labs(
        x = "% Mangrove cover",
        y = "Grass shrimp abundance",
        color = "% Land",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1), n.breaks = 4) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +  theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    ##### edge plot #####
    edge <- 3
    buf <- 250
    species_code <- "PALSP"
    run_scale = "satscale"
    # name scale habitat data file 
    habx <- paste("google2022_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "satscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    
    vis.gam(best_model)
    library(mgcv)
    
    # Calculate Cook's distance for the best_model
    cooks_distances <- cooks.distance(best_model)
    
    # Examine the Cook's distances
    summary(cooks_distances)
    
    # Plot Cook's distances
    plot(cooks_distances, type = "h", 
         main = "Cook's Distance", 
         ylab = "Cook's Distance", 
         xlab = "Observation Index")
    abline(h = 4 / length(cooks_distances), col = "red", lty = 2) # Threshold
    
    threshold <- 4 / nrow(pf.env)
    
    # Identify influential observations
    influential_obs <- which(cooks_distances > threshold)
    
    # Print influential observations
    print(influential_obs)
    weights <- rep(1, nrow(pf.env)) # ADDED DUE TO HIGH OUTLIERS at obs 1 and 49 (cooks distance of 0.4 and 0.86 respectively)
    #weights[c(1)] <- 0.1  # Assign lower weight to influential observation
    
    best_model <- gam(as.formula(formula_table[[species_code]][[run_scale]]), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML", weights = weights)
    summary(best_model)
    
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    lwr22 = smooth_table["s(land_water_ratio)",4]
    edge22 <- smooth_table["s(edge_perc)",4]
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab_edg <- create_sequence_table(pf.env, "edge_perc") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    
    # predict response values over the range of variable of interest
    predict_tab_edg$predicted_count <- predict(best_model, newdata = predict_tab_edg, type = "response")
    
    
    
    
    pf.env$true_predict <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    # Compute residuals
    pf.env$residuals <- pf.env$PALSP - pf.env$true_predict
    
    # Apply LOWESS smoothing
    loess_fit <- loess(true_predict ~ edge_perc, data = pf.env, span = 0.5)    
    pf.env$Smoothed_Pred <- predict(loess_fit)
    # Plot observed points colored by residuals and smoothed prediction curve
    for (col in predictor_columns) {
      predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'prediction.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    ggplot(pf.env, aes(x = edge_perc, y = PALSP)) +
      geom_point(aes(color = edge_l.mangrove), size = 6, alpha = 1) + 
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      
      # Map linetype inside aes() to include Smoothed_Pred in the legend
      geom_line(aes(y = Smoothed_Pred, linetype = "Smoothed \nprediction"), color = "red", size = 1.2) +
      
      # Use scale_linetype_manual to ensure it appears in the legend
      scale_linetype_manual(name = "", values = c("Smoothed \nprediction" = "solid")) +
      
      labs(x = "‰ Edge", y = "Grass shrimp abundance", col = "% Mangrove cover",
           title = sprintf("Edge distance = %s m \n Buffer radius = %s m", edge, buf)) +
      scale_x_continuous(labels = percent_format(scale = 1000, suffix = "‰")) +
      
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    dev.off()
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('satscale_', species_code,'edge.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab_edg, aes(x = edge_perc*1000, y = predicted_count, col = edge_l.mangrove)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = edge_perc*1000, y = !!ensym(species_code)), size = 4) +
      annotate("text", x = mean(pf.env$edge_perc*1000), y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(edge22, 4))) +
      labs(
        x = "‰ Edge",
        y = "Grass shrimp abundance",
        col = "% Mangrove cover",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1), n.breaks = 5) +
      scale_x_continuous(labels = label_number(suffix = "‰")) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    
    
  }
  if(run_scale == "smallscale"){
    ###### mangrove #######
    # Based on GAM outputs, select a scale to plot species abundance
    edge <- 1 
    buf <- 300
    species_code <- "PALSP"
    
    # name scale habitat data file 
    habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(top_model_formula), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    # Calculate Cook's distance for the best_model
    cooks_distances <- cooks.distance(best_model)
    
    # Examine the Cook's distances
    summary(cooks_distances)
    
    # Plot Cook's distances
    plot(cooks_distances, type = "h", 
         main = "Cook's Distance", 
         ylab = "Cook's Distance", 
         xlab = "Observation Index")
    abline(h = 0.8, col = "red", lty = 2) # Threshold
    # Residual plot
    
    
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    man22<-smooth_table["s(edge_l.mangrove)",4]
    
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab <- create_sequence_table(pf.env, "edge_l.mangrove") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    # predict response values over the range of variable of interest
    predict_tab$predicted_count <- predict(best_model, newdata = predict_tab, type = "response")
    
    # Back-transform each variable in predict_tab
    for (col in predictor_columns) {
      predict_tab[[col]] <- back_transform(predict_tab[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    # Generate predictions on the original data (pf.env)
    pf.env$predicted_count <- exp(predict(best_model, newdata = pf.env, type = "link"))
    
    library(ggplot2)
    
    # Plot observed vs. predicted values
    ggplot(data = pf.env, aes(x = .data[[species_code]], y = predicted_count)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, color = "red") +
      labs(x = "Observed", y = "Predicted") +
      theme_classic()
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('smallscale_', species_code,'mangrove.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab, aes(x = edge_l.mangrove, y = predicted_count)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = edge_l.mangrove, y = !!ensym(species_code), col = (edge_perc)*1000), size = 4) +
      annotate("text", x = mean(pf.env$edge_l.mangrove) * 0.5, y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(man22, 3))) +
      labs(
        x = "% Mangrove edge",
        y = "Grass shrimp abundance",
        color = "‰ Edge",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = label_number(suffix = "‰")) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    ##### edge plot #####
    edge <- 1 
    buf <- 100
    species_code <- "PALSP"
    
    # name scale habitat data file 
    habx <- paste("combined_", "edge", edge, "_buf", buf, sep = "")
    
    # Load habitat data
    hab <-  read.csv(file.path("./landscape_analysis", "output", "smallscale", paste(habx,".csv", sep = "")))
    hab <- hab %>% mutate('edge_perc' = (edge_man + edge_mar) / (pi * buf^2))
    
    # Merge habitat data with community composition table
    pf.env <- merge(hab, comtab[, c("site_date_key", species_code[1], "salinity", "temperature", "latitude", "longitude")], by = 'site_date_key') 
    # Standardize variables and store statistics for back-transformation
    predictor_columns <- c("edge_perc", "edge_l.mangrove", "land_water_ratio")
    standardization_stats <- data.frame(variable = predictor_columns, mean = numeric(length(predictor_columns)), sd = numeric(length(predictor_columns)))
    
    for (col in predictor_columns) {
      standardization_stats[standardization_stats$variable == col, "mean"] <- mean(pf.env[[col]], na.rm = TRUE)
      standardization_stats[standardization_stats$variable == col, "sd"] <- sd(pf.env[[col]], na.rm = TRUE)
      pf.env[[col]] <- (pf.env[[col]] - standardization_stats[standardization_stats$variable == col, "mean"]) / standardization_stats[standardization_stats$variable == col, "sd"]
    }
    
    # run the model selected above at the scale chosen and save the model
    best_model <- gam(as.formula(top_model_formula), data = pf.env, na.action = na.fail, family = tw(link = "log"), method = "ML")
    summary(best_model)
    
    # Residual plot
    residuals_gam <- residuals(best_model, type = "response")
    fitted_values <- fitted(best_model)
    
    plot(fitted_values, residuals_gam, 
         xlab = "Fitted Values", ylab = "Residuals", 
         main = "Residuals vs Fitted Values")
    abline(h = 0, col = "red")
    
    
    #extract coefficients
    smooth_table<-summary(best_model)$s.table
    man22<-smooth_table["s(edge_l.mangrove)",4]
    # Produce a datatable where the variable of interest varies over it's range while all others in the model stay the same. 
    predict_tab_edg <- create_sequence_table(pf.env, "edge_perc") # variable in the function will be sequenced, all others in the datatable will be held at the mean
    
    # predict response values over the range of variable of interest
    predict_tab_edg$predicted_count <- predict(best_model, newdata = predict_tab_edg, type = "response")
    
    for (col in predictor_columns) {
      predict_tab_edg[[col]] <- back_transform(predict_tab_edg[[col]], col, standardization_stats)
      pf.env[[col]] <- back_transform(pf.env[[col]], col, standardization_stats)
    }
    
    
    png(file.path('./species_analysis/output/plots/single_scale', paste('smallscale_', species_code,'edge.png', sep = "")), width = 8, height = 4, units = 'in', res = 300)
    
    # plot
    ggplot(predict_tab_edg, aes(x = edge_perc*1000, y = predicted_count, col = edge_l.mangrove)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(data = pf.env, aes(x = edge_perc*1000, y = !!ensym(species_code)), size = 4) +
      annotate("text", x = mean(pf.env$edge_perc*1000) * 0.5, y = max(pf.env[,as.character(species_code)]*0.8),
               label = paste("p =", round(edge22, 3))) +
      labs(
        x = "‰ Edge",
        y = "Grass shrimp abundance",
        col = "% Mangrove edge",
        title = sprintf("Edge distance = %s m \n Buffer radius = %s m",
                        edge, 
                        buf)        
      ) +
      scale_color_binned(type = 'viridis', labels = percent_format(accuracy = 1)) +
      scale_x_continuous(labels = label_number(suffix = "‰")) +
      theme_classic() +
      theme(
        text = element_text(size = 16),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      )
    
    
    dev.off()
    
    
    
  }
  
}

