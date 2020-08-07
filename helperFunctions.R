####################################################################################################################
######    MCL-1 gains occur with high frequency in lung adenocarcinoma and can be targeted therapeutically     #####
######                                          Helper Functions                                               #####
####################################################################################################################
# code written by Nicholas McGranahan (nicholas.mcgranahan.10@ucl.ac.uk) and Michelle Dietzen (m.dietzen@ucl.ac.uk)

# function to combine following segments with the same copynumber
shrink.seg <- function(chr.seg,colNameToSelect="cpnToUse"){
  new.chr <- matrix(0,0,6)
  for(i in unique(chr.seg$sample)){
    tmp <- chr.seg[chr.seg$sample %in% i,]
    if(any(duplicated(tmp[,colNameToSelect]))){
      seg.class <- c(1)
      for(j in 2:nrow(tmp)){
        if(tmp[(j-1),colNameToSelect] == tmp[j,colNameToSelect]){
          seg.class <- c(seg.class, seg.class[j-1])
        }
        if(tmp[(j-1),colNameToSelect] != tmp[j,colNameToSelect]){
          seg.class <- c(seg.class, seg.class[j-1]+1)
        }
      }
      for(j in unique(seg.class)){
        tmp[seg.class %in% j,]$endpos <- max(tmp[seg.class %in% j,]$endpos)
      }
      tmp<- tmp[!duplicated(seg.class),]
    }
    new.chr <- rbind(new.chr, tmp)
  }
  colnames(new.chr) <- colnames(chr.seg)
  return(new.chr)
}

shrink.seg.wrapper <- function(seg,colNameToSelect="cpnPloidyCorrected"){
  new.seg <- seg[1,]
  for(j in unique(seg$sample)){
    sample.seg <- seg[seg$sample %in% j,]
    new.sample.seg <- seg[1,]
    for(i in unique(sample.seg$chr)){
      sample.chrom.seg <- sample.seg[sample.seg$chr %in% i,,drop=F]
      if(nrow(sample.chrom.seg) > 1){
        sample.chrom.seg <- shrink.seg(sample.chrom.seg,colNameToSelect)
      }
      new.sample.seg <- rbind(new.sample.seg, sample.chrom.seg)
    }
    new.seg <- rbind(new.seg, new.sample.seg[-1,])
  }
  seg <- new.seg[-1,]
  return(seg)
}


#function to calculate minimum consistent  regions
min.cons.fn <- function(segments, chromosomes = 1:22, mc.cores=2,colToUse = "cpnPloidyCorrected"){
  #assumes the following columns (and order)
  # sample, chr, startpos, endpos
  require(parallel)
  segments <- as.matrix(segments)
  chrom.all<- c()
  start.all<- c()
  end.all<- c()
  
  for(i in chromosomes){
    #print(i)
    s.chrom <- subset(segments, as.numeric(segments[, 2]) == i)
    start   <- unique(as.numeric(s.chrom[, 3]))
    end     <- unique(as.numeric(s.chrom[, 4]))
    end1    <- 0
    start.final <- c()
    end.final   <- c()
    
    while(end1 != max(end)){
      if(sum(end1>as.numeric(s.chrom[,3]) & end1<as.numeric(s.chrom[,4]))>0){
        start1          <- end1+1
      } else {
        start1          <- min(start[which(start > end1)])
      }
      end1 <- min(end[end >= start1], start[start>start1])
      if(end1%in%c('Inf')){
        end1      <- min(end[end >= start1])
      }
      start.final <- c(start.final, start1)
      end.final   <- c(end.final, end1)
    }
    
    chrom.all <- c(chrom.all, rep(i, length(start.final)))
    start.all <- c(start.all, start.final)
    end.all   <- c(end.all, end.final)
  }
  
  mat.pos <- matrix(NA, length(start.all), nc = 3)
  colnames(mat.pos)<- c("Chrom", "Start", "End")
  mat.pos<- cbind(as.numeric(chrom.all), as.numeric(start.all), as.numeric(end.all))

  tmp.out <- mclapply(unique(mat.pos[,1]), function(x){
    t(apply(mat.pos[mat.pos[,1] %in% x,,drop=F], 1, fun.mapcopy.pos, seg = segments,colToUse=colToUse))}, mc.cores=mc.cores)
  mat.tmp.cn <- do.call(rbind, tmp.out)
  mat.final.cn<- cbind(mat.pos, mat.tmp.cn)
  colnames(mat.final.cn)<- c("Chr", "Start", "End", as.character(unique(segments[,1])))
  return(mat.final.cn)
}

#function to calculate frequency of clonal and subclonal copy number events
freq_clonalcn.fn <- function(min.cons.seg, sample.names, type, gain.threshold=log2(2.5/2), loss.threshold=log2(1.5/2),min.gain.threshold=gain.threshold, min.loss.threshold=loss.threshold, min.samples=2){ 
  
  samples <- sort(unique(substr(sample.names,3,8)))
  gain.min.cons <- loss.min.cons <- gain.trunk.min.cons <- gain.branch.min.cons <-loss.trunk.min.cons <-loss.branch.min.cons <- cbind(min.cons.seg[,1:3],matrix(0, nrow=nrow(min.cons.seg), ncol=length(samples)))
  colnames(loss.min.cons) <-colnames(gain.min.cons) <-colnames(gain.trunk.min.cons)<-colnames(gain.branch.min.cons) <- colnames(loss.trunk.min.cons) <-colnames(loss.branch.min.cons) <- c(colnames(gain.trunk.min.cons)[1:3], sort(samples))
  
  for(i in 4:ncol(gain.trunk.min.cons)){
    min.cons.samples <- min.cons.seg[,substr(colnames(min.cons.seg),3,8) %in% colnames(gain.trunk.min.cons)[i],drop=F]
    
    clonalgain <- t(apply(min.cons.samples,1,function(x){
      tmp <- x
      tmp[tmp>=gain.threshold] <- 1
      tmp <- tmp[!is.na(tmp)]
      if(any(tmp >= 1)){
        tmp[tmp>= min.gain.threshold] <- 1
      }
      tmp[tmp < 1] <- 0
      c(sum(tmp),length(tmp))
    }))
    clonalloss <- t(apply(min.cons.samples,1,function(x){
      tmp <- x
      tmp[tmp<= loss.threshold] <- -1
      tmp <- tmp[!is.na(tmp)]
      if(any(tmp <= -1)){
        tmp[tmp<= min.loss.threshold] <- -1
      }
      tmp[tmp > -1] <- 0
      c(abs(sum(tmp)), length(tmp))
    }))
    
    gain.branch.min.cons[clonalgain[,1] >0 & clonalgain[,1] < clonalgain[,2],i] <- 1
    loss.branch.min.cons[clonalloss[,1] >0 & clonalloss[,1] < clonalloss[,2],i] <- -1
    
    gain.trunk.min.cons[clonalgain[,1] == clonalgain[,2],i] <- 1
    loss.trunk.min.cons[clonalloss[,1] == clonalloss[,2],i] <- -1
    
    
    gain.min.cons[clonalgain[,1] >0,i] <- 1
    loss.min.cons[clonalloss[,1] >0,i] <- -1
    gain.min.cons[clonalgain[,2] < min.samples,i] <- gain.trunk.min.cons[clonalgain[,2] < min.samples,i] <- gain.branch.min.cons[clonalgain[,2] < min.samples,i] <- NA
    loss.min.cons[clonalloss[,2] < min.samples,i] <- loss.trunk.min.cons[clonalloss[,2] < min.samples,i] <-loss.branch.min.cons[clonalloss[,2] < min.samples,i] <- NA
  }
  
  if(type=='all' ){
    outlist <- list(gain.min.cons, loss.min.cons, gain.trunk.min.cons, loss.trunk.min.cons, gain.branch.min.cons, loss.branch.min.cons)
    names(outlist) <- c('all_gain','all_loss','trunk_gain','trunk_loss','branch_gain','branch_loss')
  }
  if(type=='clonal' ){
    outlist <- list(gain.trunk.min.cons, loss.trunk.min.cons)
  }
  if(type=='subclonal' ){
    outlist <- list(gain.branch.min.cons, loss.branch.min.cons)
  }
  return(outlist)
}



#function to extract gene position from position file
get.genepos <- function(gene=NULL, genepos=NULL, posfile, load.genepos=FALSE, return.genepos=FALSE){
  if(load.genepos){
    genepos <- read.table(posfile,sep='\t',header=T,as.is=T)
    genepos[,'chrom'] <- sub('chr','',genepos[,'chrom'])
    genepos[,'chrom'] <- sub('X',23,genepos[,'chrom'])
    genepos[,'chrom'] <- sub('Y',24,genepos[,'chrom'])
    genepos <- genepos[genepos[,'chrom'] %in% 1:24,]
  }
  if(return.genepos){
    return(genepos)
  }
  if(!return.genepos){
    if(is.null(genepos) | is.null(gene)){ stop('\ngene and genepos must be given\n')}
    gene_loc <- genepos[genepos[,'name2'] %in% gene,,drop=F]
    nchr <- length(unique(gene_loc[,'chrom']))	
    if(nrow(gene_loc) > 0 & nchr ==1){
      gene_chr <- as.numeric(as.character(names(sort(table(gene_loc[,'chrom'])))))
      gene_loc <- gene_loc[gene_loc[,'chrom'] %in% gene_chr,]
      gene_start <- min(gene_loc[,'txStart'])
      gene_end <- max(gene_loc[,'cdsEnd'])	
      genepos_out <- setNames(c(gene_chr, gene_start, gene_end), c('Chr','Start','End'))
    } else {
      genepos_out <- setNames(c(NA, NA, NA), c('Chr','Start','End'))
    }
    return(genepos_out)
  }
}

# function to extract copy number of specific position (e.g. gene) from copy number segment table
get.cn <- function(seg, chr, start, stop = start, mean.seg=TRUE, adjusted.seg=FALSE){
  #20140519: added adjusted.seg: if true, will calculate the mean segment size for just the region of interest, ie, it will weigh the average copy number by the size of the region affected
  
  seg <- seg[seg[,2] == chr,]
  seg <- seg[seg[,3] <= stop & seg[,4] >= start,]
  if(mean.seg==TRUE){
    dup.names <- unique(as.character(seg[,1][duplicated(seg[,1])]))
    for(i in dup.names){
      if(!adjusted.seg){
        tmp <- seg[seg[,1] == i,]
        mean.vect <- vector()
        for(j in 1:nrow(tmp)){
          mean.vect <- c(mean.vect, rep(tmp[j,6], tmp[j,5]))
        }
        tmp[1,4] <- tmp[nrow(tmp),4]
        tmp[1,5] <- sum(tmp[,5])
        tmp[1,6] <- mean(mean.vect)
      }
      if(adjusted.seg){
        tmp <- seg[seg[,1] == i,]
        tmp.start <- tmp[1,3]
        tmp.end <- tmp[nrow(tmp),4]
        tmp[1,3] <- start
        tmp[nrow(tmp),4] <- stop
        #round down to 10 kb, to increase speed (a lot!!)
        tmp[,3] <- floor(tmp[,3]/10000)
        tmp[,4] <- ceiling(tmp[,4]/10000)
        mean.vect <- vector()
        for(j in 1:nrow(tmp)){
          mean.vect <- c(mean.vect, rep(tmp[j,6], tmp[j,4]-tmp[j,3]))
        }
        tmp[1,3] <- tmp.start
        tmp[1,4] <- tmp.end
        tmp[1,5] <- sum(tmp[,5])
        tmp[1,6] <- mean(mean.vect)
      }	
      seg[seg[,1] == i,] <- tmp
    }
    seg <- seg[!duplicated(seg[,1]),]
  }
  return(seg)
}

