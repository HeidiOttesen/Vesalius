
KNN.DNA.Unique <- function(bead, df, k){
  ## Start
  #- Make sure you have the count df and bead info from the DNA_vesalius markdown..
  #- Barcodes and coordinates from bead_locations
  barcodes <- rownames(bead@coordinates)
  bead.coords <- bead@coordinates
  coords <- cbind(bead.coords$xcoord, bead.coords$ycoord)
  
  #Perform KNN:
  cat( paste(Sys.time()," Running K-Nearest Neighbor with k =",k, "number of neighbors", "\n", "Number of input barcodes:", length(unique(barcodes)), "\n"))
  coords.knearneigh <- knearneigh(coords, k = k)
  knnIx <- coords.knearneigh$nn
  

  
  ## Unique barcode knn list
  #- take the first barcode from the list and all of its neighbours (their barcode indexes) 
  #- put the string barcodes in the new list if the barcodes still exists in bc and haven't been used yet.
  #- Remove those barcodes from your barcodes vector so that each barcode is only used once
  #- Lots of groups with only 1/2/3... neighbors - not full groups
  #- There are also sometimes leftovers not placed in any group
    
  bc <- barcodes
  knn.bc1 <- vector(mode="list", length=l) #holds unique barcodes in groups of up to k
  #doubles <- vector(mode="list", length=l) #Probably don't need this..
  knn.bc2 <- vector(mode="list", length=l) #unique again but without NA elements (Still null groups)
  
  #Simpler? Try next time for unique!:
  cat(paste(Sys.time()," Using only unique barcodes in KNN groups", "\n"))
  for(i in seq_along(barcodes)){
    b <- barcodes[knnIx[i,j]]
    if(b %in% bc){
      knn.bc1[[i]] <- b
    }
    bc <- bc[!bc %in% b]
    knn.bc2[[i]] <- x[!is.na(knn.bc1[[i]])]
  }
  
  
  
  cat(paste(Sys.time()," Number of barcodes not included: ", length(bc), "\n"))
  
  
  ## Find median coordinate within knn group
  #- iterate over groups of knn 
  #- stores their corresponding x and y coordinates in vector t
  #- Stores the median x and y for each column in new list - comb
  #- makes a new list of new barcode names, bc.g - the first barcode in each group
  
  i <- 1
  l <- length(knn.bc2)
  comb <- matrix(0, nrow = l, ncol = ncol(bead.coords))
  bc.g <- vector(length=l)
  cat(paste(Sys.time()," Finding median coordinate for each knn group", "\n"))
  while(i <= l){
    tmp <- knn.bc2[[i]]
    t <- bead.coords[tmp,]
    ltmp <- length(tmp)
    comb[i,] <- apply(t, 2, median)
    if(ltmp > 0){
      bc.g[[i]] <- tmp[[1]]
    }
    i <- i + 1 
  }
  
  length(unique(bc.g)) #15737
  length(bc.g) #38313
 
  rownames(comb) <- bc.g
  colnames(comb) <- c("xcoord", "ycoord")
  
  #Remove NA values and empty groups
  
  row.has.na <- apply(comb, 1, function(x){any(is.na(x))})
  comb <- comb[!row.has.na,]
  knn.bc <- knn.bc2[!sapply(knn.bc2,is.null)]
  
  ## Summing counts over knn-groups:
  #- Iterate over knn groups 
  #- take the counts for those barcodes for all bins from the sparsematrix
  #- df - Sparse matrix vector from the DNA vesalius script - Barcodes as columns, bins as rows.
  #- Sum the counts for each group and each bin 
  
  spMtrx <- df
  i <- 1
  
  l <- length(knn.bc)
  cat(paste(Sys.time()," Number of KNN groups with only unique barcodes: ", l, "\n"))
  grMtrx <- matrix(0, nrow = nrow(spMtrx), ncol = l)
  cat( paste(Sys.time()," Summing the counts for each group and bin", "\n"))
  while(i <= l){
    tmp <- knn.bc[[i]]
    t <- spMtrx[,tmp]
    ltmp <- length(tmp)
    if(ltmp > 1){
      grMtrx[,i] <- rowSums(t)
    }else if(ltmp == 1){
      grMtrx[,i] <- t
    }
    i <- i + 1 
  }
  
  #Naming rows using new list of barcode names (first barcode of each knn group)
  bc.comb <- rownames(comb)
  

  # If reading slideseq bead info from file:
  #fname <- paste0(in.path, "knn/", alias, "_knn", k)
  #write.csv(comb, file = paste0(fname, ".bead_locations.csv"), row.names = TRUE, quote = FALSE)
  #cat(paste(Sys.time()," New bead location file saved as ", paste0(fname, ".bead_locations.csv"), "\n"))
  #ss <- ReadSlideSeq(paste0(fname, ".bead_locations.csv"))
  
  # If creating slideseq bead info as object in list instead:
  comb <- as.data.frame(comb)
  ss <- new(Class = 'SlideSeq',assay = "Spatial",
            coordinates = comb[,c("xcoord","ycoord")])
  rownames(ss@coordinates) <- rownames(comb)
  
  colnames(grMtrx) <- bc.comb
  rownames(grMtrx) <- bins
  knnSpMtx <- Matrix(grMtrx, sparse = TRUE)
  return(list(knnSpMtx, ss))
}








KNN.DNA <- function(bead, df, k){
  # Input:  Count matrix as dataframe 
  #         Bead info from the DNA_vesalius markdown..
  # Barcodes and coordinates from bead_locations
  barcodes <- rownames(bead@coordinates)
  bead.coords <- bead@coordinates
  coords <- cbind(bead.coords$xcoord, bead.coords$ycoord)
  
  #Perform KNN:
  cat(paste(Sys.time()," Running K-Nearest Neighbor with k =",k, "number of neighbors \n on", length(barcodes), "barcoded bead locations", "\n"))
  coords.knearneigh <- knearneigh(coords, k = k)
  knnIx <- coords.knearneigh$nn
  
  #Replace barcode index number with barcode string
  knn.bc <- vector(mode="list", length=length(barcodes)) #holds unique barcodes in groups of up to k
  cat(paste(Sys.time()," Change to barcode names (instead of barcode index numbers)", "\n"))
  for(i in seq_along(barcodes)){
    knn.bc[[i]] <- barcodes[knnIx[i,]]
  }
  
  
  ## Summing counts over knn-groups:
  #- Iterate over knn groups - take their barcodes 
  #- take the counts for those barcodes for all bins from the sparsematrix df
  #- df - Sparse matrix vector from the DNA_vesalius markdown script - Barcodes as columns, bins as rows.
  #- Sum the counts for each group and each bin 
  grMtrx <- matrix(0, nrow = nrow(df), ncol = ncol(df))
  cat( paste(Sys.time()," Summing the counts for each KNN group and genomic bin", "\n"))
  for(i in seq_along(knn.bc)){
    grMtrx[,i] <- rowSums(df[,knn.bc[[i]]])
  }
  
  # If creating bead info and count matrix as object in list instead of exporting to files:
  colnames(grMtrx) <- rownames(bead.coords)
  rownames(grMtrx) <- bins
  #knnSpMtx <- Matrix(grMtrx, sparse = TRUE)
  return(list(grMtrx, bead.coords))
}