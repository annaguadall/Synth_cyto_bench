## This script defines and exports functions as R objects

# Min-max normalization
normalize <- function(x) { return ((x - min(x)) / (max(x) - min(x))) }

# Matching labels and predictions according the F1 score and merging populations that have been clustered together
matching <- function(clustering, labels, threshold){
    c <- table(clustering, labels) # Contingency matrix
    r <- prop.table(c, 2)          # Recall matrix
    p <- prop.table(c, 1)          # Precision matrix
    f <- 2 * p * r / (p + r)       # F1 matrix
    f[is.na(f)] <-  0              # replacing NAs by zeros
    
    # Finding the cell population (columns) 
    # with a higher score F1 for each cluster (rows):
    m <- apply(f, 1, function(x){which.max(x)})
    m <- as.list(m)
    
    # Normalize column maximums
    col_norm <- apply(f, 2, normalize)
    
    # Col-norm-F1 column equivalent maximums
    col_m <- apply(col_norm, 2,  function(x){ which(max(x)  - x < threshold & x != 0) }  )
    col_m <- as.list(col_m)
    
    # Normalize row values (row-norm-F1)
    row_norm <- apply(f, 1, normalize) # why is the result transposed?
    row_norm <- t(row_norm)
    
    # Row-norm-F1 row equivalent maximums
    row_m <- apply( row_norm, 1, function(x){ which(max(x) - x < threshold & x != 0) }  )
    row_m <- as.list(row_m)
    
    # MERGING POPULATIONS (IN CASE ONE CLUSTER MATCHES MANY POPULATIONS)
    merged_labels <- as.character(labels)
    
    count_merged_pops <- rep(1, length(labels))
    
    for (i in 1:length(labels)){
        for (j in unlist(col_m)){ 
            if(labels[i] %in% unlist(names(row_m[[j]]))){
                if(unlist(names(row_m[[j]]))[1] == names(col_m[unlist(names(row_m[[j]]))[1]])){
                    merged_labels[i] <- paste( unlist( names(row_m[[j]]) ), collapse = "-" )
                    count_merged_pops[i] <- length(unlist( names(row_m[[j]]) ))
                }
            }
        }
    }
    
    
    # MATCHING CLUSTERS TO LABELS
    # Empty list 
    matched_merged_clusters <- rep("NA", length(clustering))
    matched_clusters <- rep("NA", length(clustering))
    
    # Replacing the numbers of the clusters by the names of the cell types:
    for(i in 1:length(clustering)){
        cat(i, " clusters,")
        for(j in 1:length(row_m)){ 
            cat(j, "row m\n")
            if(clustering[i] == names(row_m)[j]){ # if  a cluster number equals a row maximum...
                # ... give it the name of the absolute row maximum
                matched_clusters[i] <- levels(labels)[m[[j]]]
                # ... give it the name of the corresponding populations,
                # merging the name to the other equivalent row maximums, if there are
                matched_merged_clusters[i] <- paste(levels(labels)[as.numeric(unlist(row_m[[j]]))], 
                                                    collapse = "-")
            }
        }
    }
    
    # Number of PARTTIONS (matched clusters) 
    partitions <- length(table(matched_clusters))
    
    # Factorize matched (merged) clusters and labels 
    merged_labels <- as.factor(merged_labels)
    matched_merged_clusters <- factor(matched_merged_clusters, levels = levels(merged_labels)) 
    matched_clusters <- factor(matched_clusters, levels = levels(labels)) 
    # same level order than merged_labels
    # adds labels that have not been predicted
    
    matched <- list("c" = c, "f" = f, "m" = m, "merged_labels" = merged_labels,
                    "clusters" = matched_clusters, "merged_clusters" = matched_merged_clusters,
                    "partitions" = partitions, "count_merged_pops" = count_merged_pops)
    return(matched)
}

# Computing DIFFERENT MEAN F1 scores
mean_f1 <- function(cm, cm_merged, merged_labels, count_merged_pops){
    # NOT MERGED, ALL
    f1_list <- cm$byClass[,"F1"]
    f1_zeros <- f1_list
    f1_zeros <- ifelse(is.na(f1_list), 0, f1_list)
    mf1 <- mean(f1_zeros) # MEAN
    
    # NOT MERGED, ALL, WEIGHTED
    pop_freq <- prop.table(table(labels))
    mf1_w <- t(pop_freq) %*% f1_zeros 
    
    # NOT MERGED, ALL, INVERSED WEIGHTS
    inv_freq <- prop.table(1/table(labels))
    mf1_i_w <- t(inv_freq) %*% f1_zeros 
    
    # MERGED, ALL
    merged_f1_list <- cm_merged$byClass[,"F1"]
    merged_f1_zeros <- merged_f1_list
    merged_f1_zeros <- ifelse(is.na(merged_f1_list), 0, merged_f1_list)
    mf1_m <- mean(merged_f1_zeros) # MEAN
    
    # MERGED, CORRECTED
    merged_table <- table(merged_labels, count_merged_pops)
    merged_table <- cbind(merged_table, correction = 0)
    
    for(i in 1:nrow(merged_table)){
        for(j in 1:(ncol(merged_table)-1)){
            if(merged_table[i, j] != 0){
                merged_table[i, "correction"] <- merged_f1_zeros[i]/as.numeric(colnames(merged_table)[j])
            }
        }
    }
    
    mf1_m_all_corrected <- mean(merged_table[,"correction"]) # CORRECTED MEAN
    
    f1 <- list("mf1" = mf1, "mf1_w" = mf1_w, "mf1_i_w" = mf1_i_w,
               "mf1_m" = mf1_m, "mf1_m_c" = mf1_m_all_corrected)
    
    return(f1)
}
