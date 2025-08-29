#' Ternary diagrams with the CI function results
#'
#' Display the CI prediction for the selected number 
#'
#' @param data Data frame containing the results from the CI_Method function
#' @param tracers Number of tracers to display
#' @param Means Boolean indicating whether to use mean and SD data (TRUE) or raw data (FALSE)
#' @param n_col Number of columns in which to display the tracers
#' @param n_row Number of rows in which to display the tracers
#' @param cluster Boolean indicating whether to order the triangles based on a clustering function
#' @param dendro_lvls Number of branch levels in the clustering dendrogram
#' @param vrtl_sol Boolean indicating whether to plot a theoretical solution
#' @param vrtl Numeric vector of values for the theoretical solution
#' 
#' @return Ternary diagrams from the individual tracer predictions
#'
#' @export
#'
Ternary_diagram <- function(data, tracers = 1:2, n_row = 1, n_col = 2, Means = TRUE, cluster = TRUE, 
                            dendro_lvls = 5, vrtl_sol = FALSE, vrtl = c(0.05, 0.9, 0.05)) {
  
  sources <- nrow(inputSource(data[[length(data) - 1]]))
 
  if (sources == 3) {
    if (Means == TRUE) {
      if (cluster == TRUE){
      ######### Cluster ########################################################  
      df <- as.data.frame(data[(length(data)-1)])
      colnames(df) <- sub("^DB\\.", "", colnames(df))

      df <- df[c(1:nrow(df)-1),]
      
      # Generate synthetic data for each element
      synthetic_data <- df %>%
        rowwise() %>%
        mutate(
          across(starts_with("D"), ~ list(rnorm(n[1], mean = get(sub("^D", "", cur_column())), sd = .)),
                 .names = "{.col}_synthetic")
        ) %>%
        ungroup() %>%
        # Select synthetic columns
        select(ends_with("_synthetic")) %>%
        # Unnest synthetic columns
        unnest(cols = everything()) %>%
        # Rename columns to remove "D" prefix and "_synthetic" suffix
        rename_with(~ gsub("^D|_synthetic$", "", .x))  

      # Select numeric columns (excluding the first three columns)
      numeric_columns <- df %>% select(-(1:2)) %>% select_if(is.numeric)
      
      # Compute the correlation matrix
      correlation_matrix <- cor(synthetic_data, use = "pairwise.complete.obs")
      
      # Create a hierarchical clustering dendrogram
      hc <- hclust(as.dist(1 - correlation_matrix), method = "ward.D2")#"complete"
      
      # Cut the dendrogram at different levels
      clusters_levels <- cutree(hc, k = dendro_lvls) 
      
      # Set the color palette for the correlation plot
      color_palette <- colorRampPalette(c("darkgreen", "orange", "purple"))(dendro_lvls)
      
      data[1:(length(data)-2)] <- data[hc$order]
      
      # Create a logical vector for column selection
      cols_to_keep <- !(seq_len(ncol(data[[length(data) - 1]])) %in% c(1, 2, ncol(data[[length(data) - 1]])) | grepl("^D", colnames(data[[length(data) - 1]])))
      
      # Get the selected column names
      t_names <- colnames(data[[length(data) - 1]])[cols_to_keep]
      t_names[1:length(t_names)] <- t_names[hc$order]
      
      } else{
        t_names <- colnames(data[[length(data) - 1]][3:(ncol(data[[length(data) - 1]]) - 3)])
        # Set a default color (e.g., blue) for non-clustered plots
        color_palette <- rep("blue", length(tracers))
      }  
      ################################################################################  

    } else { # Means == F
        if (cluster == TRUE){
      ######### Cluster ########################################################  
      df <- as.data.frame(data[(length(data)-1)])
      colnames(df) <- sub("^DB\\.", "", colnames(df))
      
      # Select numeric columns (excluding the first three columns)
      numeric_columns <- df %>% select(-(1:2)) %>% select_if(is.numeric)
      
      # Compute the correlation matrix
      correlation_matrix <- cor(numeric_columns, use = "pairwise.complete.obs")
      
      # Create a hierarchical clustering dendrogram
      hc <- hclust(as.dist(1 - correlation_matrix), method = "ward.D2")#"complete"
      
      # Cut the dendrogram at different levels
      clusters_levels <- cutree(hc, k = dendro_lvls) 
      
      # Set the color palette for the correlation plot
      color_palette <- colorRampPalette(c("darkgreen", "orange", "purple"))(dendro_lvls)
      
      data[1:(length(data)-2)] <- data[hc$order]
      
      t_names <- colnames(data[[length(data) - 1]][3:ncol(data[[length(data) - 1]])])
      t_names[1:length(t_names)] <- t_names[hc$order]
      
        } else {
          t_names <- colnames(data[[length(data) - 1]][3:ncol(data[[length(data) - 1]])])
          # Set a default color (e.g., blue) for non-clustered plots
          color_palette <- rep("blue", length(tracers))
        }
      ################################################################################  
    } 
    
    plots <- list()
    par(mar = c(0, 0.1, 3, 0.1), mfcol = c(n_row, n_col))
    for (i2 in tracers) {
      result <- data[[i2]]
        x <- result$w.S1
        y <- result$w.S2
        z <- result$w.S3
        
        x_vrt <- vrtl[1]
        y_vrt <- vrtl[2]
        z_vrt <- vrtl[3]
        
        test <- as.data.frame(cbind(x, y, z))
        test_vrt <- as.data.frame(cbind(x_vrt, y_vrt, z_vrt))
        
        plots[[i2]] <- TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
        if (cluster == T) {
            TernaryPoints(test, col = alpha(color_palette[clusters_levels[t_names[[i2]]]], 0.5), cex = 0.1)
        } else {
            TernaryPoints(test, col = alpha('blue', 0.5), cex = 0.1)
        }
           if (vrtl_sol) {TernaryPoints(test_vrt, col = alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
           }
         title(t_names[i2], cex.main = 3.5)
      }
  } else if (sources == 4) {
      if (Means == TRUE) {
        if (cluster == TRUE) {
        ######### Cluster ########################################################  
        df <- as.data.frame(data[(length(data)-1)])
        colnames(df) <- sub("^DB\\.", "", colnames(df))
        
        df <- df[c(1:nrow(df)-1),]
        
        # Generate synthetic data for each element
        synthetic_data <- df %>%
          rowwise() %>%
          mutate(
            across(starts_with("D"), ~ list(rnorm(n[1], mean = get(sub("^D", "", cur_column())), sd = .)),
                   .names = "{.col}_synthetic")
          ) %>%
          ungroup() %>%
          # Select synthetic columns
          select(ends_with("_synthetic")) %>%
          # Unnest synthetic columns
          unnest(cols = everything()) %>%
          # Rename columns to remove "D" prefix and "_synthetic" suffix
          rename_with(~ gsub("^D|_synthetic$", "", .x))  
        
        # Select numeric columns (excluding the first three columns)
        numeric_columns <- df %>% select(-(1:2)) %>% select_if(is.numeric)
        
        # Compute the correlation matrix
        correlation_matrix <- cor(synthetic_data, use = "pairwise.complete.obs")
        
        # Create a hierarchical clustering dendrogram
        hc <- hclust(as.dist(1 - correlation_matrix), method = "ward.D2")#"complete"
        
        # Cut the dendrogram at different levels
        clusters_levels <- cutree(hc, k = dendro_lvls) 
        
        # Set the color palette for the correlation plot
        color_palette <- colorRampPalette(c("darkgreen", "orange", "purple"))(dendro_lvls)
        
        data[1:(length(data)-2)] <- data[hc$order]
        
        ################################################################################  
        # Create a logical vector for column selection
        cols_to_keep <- !(seq_len(ncol(data[[length(data) - 1]])) %in% c(1, 2, ncol(data[[length(data) - 1]])) | grepl("^D", colnames(data[[length(data) - 1]])))
        
        # Get the selected column names
        t_names <- colnames(data[[length(data) - 1]])[cols_to_keep]
        t_names[1:length(t_names)] <- t_names[hc$order]
        
        } else{ # cluster else
        t_names <- colnames(data[[length(data) - 1]][3:(ncol(data[[length(data) - 1]]) - 3)])
        # Set a default color (e.g., blue) for non-clustered plots
        color_palette <- rep("blue", length(tracers))
       }                    
      } else { # Means else
        
         if (cluster == TRUE){
          ######### Cluster ########################################################  
          df <- as.data.frame(data[(length(data)-1)])
          colnames(df) <- sub("^DB\\.", "", colnames(df))
          
          # Select numeric columns (excluding the first three columns)
          numeric_columns <- df %>% select(-(1:2)) %>% select_if(is.numeric)
          
          # Compute the correlation matrix
          correlation_matrix <- cor(numeric_columns, use = "pairwise.complete.obs")
          
          # Create a hierarchical clustering dendrogram
          hc <- hclust(as.dist(1 - correlation_matrix), method = "ward.D2")#"complete"
          
          # Cut the dendrogram at different levels
          clusters_levels <- cutree(hc, k = dendro_lvls) 
          
          # Set the color palette for the correlation plot
          color_palette <- colorRampPalette(c("darkgreen", "orange", "purple"))(dendro_lvls)
          
          data[1:(length(data)-2)] <- data[hc$order]
          
          ################################################################################  
          
          t_names <- colnames(data[[length(data) - 1]][3:ncol(data[[length(data) - 1]])])
          t_names[1:length(t_names)] <- t_names[hc$order]
        } else {
          t_names <- colnames(data[[length(data) - 1]][3:ncol(data[[length(data) - 1]])])
          # Set a default color (e.g., blue) for non-clustered plots
          color_palette <- rep("blue", length(tracers))
        }
      }
        
    plots <- list()
    par(mar = c(0, 0.1, 3, 0.1), mfcol = c(6, n_col))
    for (i2 in tracers) {
      result <- data[[i2]]
      
      # Original triangles
      x <- result$w.S1
      y <- result$w.S2
      z <- result$w.S3 + result$w.S4
      
      x1 <- result$w.S2
      y1 <- result$w.S3
      z1 <- result$w.S1 + result$w.S4
      
      x2 <- result$w.S3
      y2 <- result$w.S4
      z2 <- result$w.S1 + result$w.S2
      
      x3 <- result$w.S4
      y3 <- result$w.S1
      z3 <- result$w.S2 + result$w.S3
      
      x4 <- result$w.S1
      y4 <- result$w.S3
      z4 <- result$w.S2 + result$w.S4 
      
      x5 <- result$w.S2
      y5 <- result$w.S4
      z5 <- result$w.S1 + result$w.S3
      
      test <- as.data.frame(cbind(x, y, z))
      test1 <- as.data.frame(cbind(x1, y1, z1))
      test2 <- as.data.frame(cbind(x2, y2, z2))
      test3 <- as.data.frame(cbind(x3, y3, z3))
      test4 <- as.data.frame(cbind(x4, y4, z4))
      test5 <- as.data.frame(cbind(x5, y5, z5))
      
      # Virtual sample
      x_vrt <- vrtl[1]
      y_vrt <- vrtl[2]
      z_vrt <- vrtl[3] + vrtl[4]
      test_vrt <- as.data.frame(cbind(x_vrt, y_vrt, z_vrt))
      x_vrt1 <- vrtl[2]
      y_vrt1 <- vrtl[3]
      z_vrt1 <- vrtl[1] + vrtl[4]
      test_vrt1 <- as.data.frame(cbind(x_vrt1, y_vrt1, z_vrt1))
      x_vrt2 <- vrtl[3]
      y_vrt2 <- vrtl[4]
      z_vrt2 <- vrtl[1] + vrtl[2]
      test_vrt2 <- as.data.frame(cbind(x_vrt2, y_vrt2, z_vrt2))
      x_vrt3 <- vrtl[4]
      y_vrt3 <- vrtl[1]
      z_vrt3 <- vrtl[2] + vrtl[3]
      test_vrt3 <- as.data.frame(cbind(x_vrt3, y_vrt3, z_vrt3))
      x_vrt4 <- vrtl[1]
      y_vrt4 <- vrtl[3]
      z_vrt4 <- vrtl[2] + vrtl[4]
      test_vrt4 <- as.data.frame(cbind(x_vrt4, y_vrt4, z_vrt4))
      x_vrt5 <- vrtl[2]
      y_vrt5 <- vrtl[4]
      z_vrt5 <- vrtl[1] + vrtl[3]
      test_vrt5 <- as.data.frame(cbind(x_vrt5, y_vrt5, z_vrt5))
      
      # Plotting
      # par(mfcol=c(6, n_col))  # Set up multi-plot layout
      # Plot original triangles
      plots[[i2]] <- TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      if (cluster == T) {
        TernaryPoints(test, col = if (cluster) alpha(color_palette[clusters_levels[t_names[[i2]]]], 0.5) else alpha(color_palette[i2], 0.5), cex = 0.1)
      } else {
        TernaryPoints(test, col = alpha('blue', 0.5), cex = 0.1)
      }
      if (vrtl_sol) {
        TernaryPoints(test_vrt,  col = alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
      title(t_names[i2], cex.main = 3.5)
      
      TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      if (cluster == T) {
        TernaryPoints(test1, col = if (cluster) alpha(color_palette[clusters_levels[t_names[[i2]]]], 0.5) else alpha(color_palette[i2], 0.5), cex = 0.1)
      } else {
        TernaryPoints(test1, col = alpha('blue', 0.5), cex = 0.1)
      }
      if (vrtl_sol) {
        TernaryPoints(test_vrt1,  col = alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
      TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      if (cluster == T) {
        TernaryPoints(test2, col = if (cluster) alpha(color_palette[clusters_levels[t_names[[i2]]]], 0.5) else alpha(color_palette[i2], 0.5), cex = 0.1)
      } else {
        TernaryPoints(test2, col = alpha('blue', 0.5), cex = 0.1)
      }
      if (vrtl_sol) {
        TernaryPoints(test_vrt2,  col = alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
      TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      if (cluster == T) {
        TernaryPoints(test3, col = if (cluster) alpha(color_palette[clusters_levels[t_names[[i2]]]], 0.5) else alpha(color_palette[i2], 0.5), cex = 0.1)
      } else {
        TernaryPoints(test3, col = alpha('blue', 0.5), cex = 0.1)
      }
      if (vrtl_sol) {
        TernaryPoints(test_vrt3,  col = alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
       TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      # TernaryPoints(test4, col = alpha(color_palette[clusters_levels[t_names[[i2]]]], 0.5), cex = 0.1)
       if (cluster == T) {
         TernaryPoints(test4, col = if (cluster) alpha(color_palette[clusters_levels[t_names[[i2]]]], 0.5) else alpha(color_palette[i2], 0.5), cex = 0.1)
       } else {
         TernaryPoints(test4, col = alpha('blue', 0.5), cex = 0.1)
       }
       if (vrtl_sol) {
        TernaryPoints(test_vrt4,  col = alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
      TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      if (cluster == T) {
        TernaryPoints(test5, col = if (cluster) alpha(color_palette[clusters_levels[t_names[[i2]]]], 0.5) else alpha(color_palette[i2], 0.5), cex = 0.1)
      } else {
        TernaryPoints(test5, col = alpha('blue', 0.5), cex = 0.1)
      }
      if (vrtl_sol) {
        TernaryPoints(test_vrt5,  col = alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
    }
  }
}
