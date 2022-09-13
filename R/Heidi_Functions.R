cimg.plot <- function(image, step) {
  ss <- subset(image, select=c(x, y, cc, value))
  g <- as.cimg(ss) %>% plot(main=c(step, "Array (as.cimg)"))
}


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
  bc <- barcodes    # copy of barcode list to remove from as we iterate - to only get unique
  knn.bc1 <- vector(mode="list", length=length(bc)) #holds unique barcodes in groups of up to k
  knn.bc2 <- vector(mode="list", length=length(bc)) #unique again but without NA elements (Still null groups)
  
  #Simpler? Try next time for unique!:
  cat(paste(Sys.time()," Using only unique barcodes in KNN groups", "\n"))
  for(i in seq_along(barcodes)){
    b <- barcodes[knnIx[i,]]
    for(j in seq_along(b)){
      br <- b[[j]]
      if(br %in% bc){       ## if current barcode still exists in the copy of barcodes we can add it to the new knn list
        knn.bc1[[i]][j] <- br
      }
    }
    bc <- bc[!bc %in% b] #remove the barcodes from current barcode group - keep only unused barcodes
    x <- knn.bc1[[i]]
    knn.bc2[[i]] <- x[!is.na(x)] #Remove all NA elements
  }
  
  cat(paste(Sys.time()," Number of barcodes not included: ", length(bc), "\n"))
  
  ## Find median coordinate within knn group
  #- iterate over groups of knn 
  #- stores their corresponding x and y coordinates in vector t
  #- Stores the median x and y for each column in new list - comb
  #- makes a new list of new barcode names, bc.g - the first barcode in each group
  l <- length(knn.bc2)
  comb <- matrix(0, nrow = l, ncol = ncol(bead.coords))
  bc.g <- vector(length=l)
  cat(paste(Sys.time()," Finding median coordinate for each knn group", "\n"))
  for(i in seq(l)){
    tmp <- knn.bc2[[i]]
    t <- bead.coords[tmp,]
    ltmp <- length(tmp)
    comb[i,] <- apply(t, 2, median)
    if(ltmp > 0){
      bc.g[[i]] <- tmp[[1]]
    }
    i <- i + 1 
  }
  
  length(unique(bc.g)) 
  length(bc.g) 
 
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
  for(i in seq_along(knn.bc)){
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
  ss <- ss@coordinates
  
  colnames(grMtrx) <- bc.comb
  rownames(grMtrx) <- bins
  knnSpMtx <- Matrix(grMtrx, sparse = TRUE)
  return(list(knnSpMtx, ss))
}






#' KNN.DNA - Perform K nearest neighbor on count matrix with coordinates
#' @param bead SlideSeq object containing barcode coordinate information
#' @param df Matrix containing the sparse counts of bins
#' @param k Number of neighbors to sum counts of

KNN.DNA <- function(bead, df, k){
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





#' chromosomeCounts - Sum counts of bins over full chromosomes
#' @param df count matrix
#' @param bin.all bin index + chromosome + chr location information

## Replace bin index with chromosome number - nope they are unique rownames in the df matrix
# Make groups of bin indexes based on chromosome index
# Sum counts of each chromosome


chromosomeCounts <- function(df, bin.all){
  ## Summing counts across chromosomes instead of bins (per barcode)
  # bin.group a list of all the indexes within each chromosome
  # chrMatrix - the summed count matrix
  chr <- (length(unique(bin.all$chr_ind)) -1)    #Number of chromosomes minus the mitochondrial
  bin.group <- list()
  chrMatrix <- matrix(0, nrow = chr, ncol = ncol(df))
  for(i in seq(chr)){
    g <- subset(bin.all, bin.all$chr_ind==i)
    bin.group[[i]] <- g$bin_ind
    chrMatrix[i,] <- colSums(df[bin.group[[i]],])
  }
  
  # If creating bead info and count matrix as object in list instead of exporting to files:
  colnames(chrMatrix) <- colnames(df)
  rownames(chrMatrix) <- 1:chr
  return(chrMatrix) 
}







#' binLocations - Retrieve chromosome location information from bin-index number
#' @param bin.chr The top interesting bins with chromosomal location information
#' @param ref.bins Reference bins - either mouse (mm10) or human (hg19)

binLocations <- function(int.bins, ref.bins){
  chr <- data_frame()
  for(i in seq_along(int.bins$genes)){
    x <- as.list(subset(ref.bins, bin_ind==int.bins$genes[[i]]))
    chr <- bind_rows(chr, x)
  }
  chr <- chr[c("bin_start", "bin_end", "chr_ind", "bin_ind")]
  
  names(int.bins)[names(int.bins) == 'genes'] <- 'bin_ind'
  bin.chr <- merge(chr, int.bins)
  bin.chr <- bin.chr[order(bin.chr$bin_start),]
  bin.chr <- bin.chr[order(bin.chr$chr_ind),]
  return(bin.chr)
}


genesBioMart <- function(bin.chr, alias){
  library(biomaRt)
  if(grepl("human",alias)){
    mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  } else if(grepl("mouse",alias)){
    mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="nov2020.archive.ensembl.org", path="/biomart/martservice", dataset="mmusculus_gene_ensembl")
  }
  
  r <- vector()
  genes <- data.frame()
  for(i in seq_along(bin.chr$bin_ind)){ #OBS! which column/$group can i seq_along - different names in chromosome vs bin?
    r <- list(bin.chr$chr_ind[i], bin.chr$bin_start[i], bin.chr$bin_end[i])
    # Skip "go_id" to get fewer results...
    g <- getBM(c("external_gene_name", "description", "ensembl_gene_id", "gene_biotype", "chromosome_name", 
                 "start_position","end_position", "strand"),
               filters=c("chromosome_name","start","end"),
               values=r, mart=mart)
    genes <- rbind(genes, g)
  }
  genes <- genes[order(genes$start_position),]
  genes <- genes[order(genes$chromosome_name),]
  return(genes)
}


