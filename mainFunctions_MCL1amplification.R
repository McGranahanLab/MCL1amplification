####################################################################################################################
######    MCL-1 gains occur with high frequency in lung adenocarcinoma and can be targeted therapeutically     #####
######                                           Main Functions                                                #####
####################################################################################################################
# code written by Nicholas McGranahan (nicholas.mcgranahan.10@ucl.ac.uk) and Michelle Dietzen (m.dietzen@ucl.ac.uk)

source('helperFunctions.R')

##### function to calculate and plot copy number frenquencies around gene of interest #####
#--> frequency plot for whole genome was created in the same way, except for excluding the part where we subset for region of interest
freq_plot_gene.fn <- function(seg
                              , purity.threshold=0
                              , gain.threshold=log2(2.5/2)
                              , loss.threshold=log2(1.5/2)
                              , type='all'
                              , pdf_name='frequency_plot'
                              , min.gain.threshold=gain.threshold
                              , min.loss.threshold=loss.threshold
                              , min.samples=1
                              , ploidy='mean'
                              , use.raw=TRUE
                              , add_pval=TRUE
                              , GeneOfInterest ='MCL1'
                              , main=""
                              , basePairsEitherSide=50000000
                              , add.gene.pos=NULL
                              , UCSCgenes.file = 'genes.UCSC.hg19.grch37.20141112.txt' 
                              , cytobands.file = 'cytoBand.hg19.txt'){
  
  #only consider samples with specific purities
  tmp <- seg[seg[,'ACF'] >= purity.threshold & seg[,'ACF'] < 1,]
  
  #different ways to consider copy number
  if(use.raw){
    tmp[,6] <- rowSums(tmp[,c('nAraw','nBraw')])
  }
  if(ploidy=='main'){
    samp.ploidy<- calc.ploidy(tmp)
    tmp[,6] <- log2(tmp[,6]/samp.ploidy[match(tmp[,1],rownames(samp.ploidy)),2])
  }
  if(ploidy=='mean'){
    tmp[,6] <- log2(tmp[,6]/tmp[,9])
  }    
 
  #only use autosomes 
  tmp <- tmp[tmp[,2] %in% c(1:22),]
  tmp <- shrink.seg.wrapper(tmp,colNameToSelect = 'cnTotal')
  sample.names_tmp <- as.character(unique(tmp[,1]))
  
  #calculate minimum consistent regions
  min.cons.seg     <- min.cons.fn(tmp, chromosomes=1:22, colToUse = 'cnTotal') 
  min.cons.seg[min.cons.seg == '-Inf'] <- -3
  
  #calculate frequency of clonal and subclonal copy number events
  min.cons.out <- freq_clonalcn.fn(min.cons.seg, sample.names_tmp, type, gain.threshold, loss.threshold, min.gain.threshold,min.loss.threshold, min.samples=min.samples)
  rownames(min.cons.out[[2]]) <- rownames(min.cons.out[[1]]) <- 1:nrow(min.cons.out[[1]])
  
  min.cons.seg.gain <- min.cons.out[[1]]
  min.cons.seg.loss <- min.cons.out[[2]]
  
  na_rms <- !apply(min.cons.seg.gain,2, function(x){all(is.na(x))}) & !apply(min.cons.seg.loss,2, function(x){all(is.na(x))})
  
  min.cons.seg.gain <- min.cons.seg.gain[,na_rms]
  min.cons.seg.loss <- min.cons.seg.loss[,na_rms]
  
  sample.names <- colnames(min.cons.seg.gain)[-c(1:3)]
  
  # calculate  significance thresholds for frequency of gain and loss using bootstrapping
  totalsum <- sum(na.omit(min.cons.seg.loss[,3] - min.cons.seg.loss[,2]))
  gain.loss.perc <- t(sapply(colnames(min.cons.seg.gain)[-c(1:3)],function(x){    
    gains <- min.cons.seg.gain[,colnames(min.cons.seg.gain) %in% x]
    sumgains <- sum(na.omit(min.cons.seg.gain[gains >0,3]- min.cons.seg.gain[gains >0,2]))
    loss <- min.cons.seg.loss[,colnames(min.cons.seg.loss) %in% x]
    sumloss <- sum(na.omit(min.cons.seg.loss[loss < 0,3]- min.cons.seg.loss[loss <0,2]))
    c(sumloss/totalsum,sumgains/totalsum)
  }))
  colnames(gain.loss.perc) <- c('Perc.loss','Perc.gain')

  fun.sample <- function(p) {return(sample(c(0, 1), size = 1,  prob = c(1-p, p)))  }
  
  perm.loss.samples <- c()  
  for(i in 1:10000)  { 
    perm.loss.samples <- c(perm.loss.samples , sum(as.numeric(sapply(gain.loss.perc[, 1], fun.sample)))  )
  }  
  perm.gain.samples <- c()  
  for(i in 1:10000)  { 
    perm.gain.samples <- c(perm.gain.samples , sum(as.numeric(sapply(gain.loss.perc[, 2], fun.sample)))  )
  }  
  
  sig.thresh.copy         <- c(quantile(perm.loss.samples, c((1-0.05)),na.rm =T),quantile(perm.gain.samples, c((1-0.05)),na.rm=T))
  sig.thresh.copy         <- rbind(sig.thresh.copy)/length(sample.names)
  
  freq.min.copy.loss     <-  apply(min.cons.seg.loss[,sample.names],1,sum,na.rm=T)
  perc.min.copy.loss     <-  freq.min.copy.loss/ncol(min.cons.seg.loss[,sample.names])
  
  #calculate p.value for each gain and loss 
  sig.test <- function (x, perm.samples) {
    p.val <- sum(perm.samples>=x) / length(perm.samples)
    p.val <- max(p.val,1e-04)
    return(p.val)
  }
  p.val.loss             <- sapply(abs(freq.min.copy.loss),sig.test,perm.loss.samples)
  freq.min.copy.gain     <-  apply(min.cons.seg.gain[,sample.names],1,sum,na.rm=T)
  perc.min.copy.gain     <-  freq.min.copy.gain/ncol(min.cons.seg.gain[,sample.names])
  p.val.gain             <-  sapply(abs(freq.min.copy.gain),sig.test,perm.gain.samples)
  
  #subset for values around gene of interest
  genepos <- get.genepos(load.genepos=TRUE, return.genepos=TRUE,posfile = UCSCgenes.file)
  a       <- get.genepos(GeneOfInterest,genepos)
  
  geneMinConsOutIndex <- which(min.cons.seg[,1] %in% a[1] & min.cons.seg[,2]<= a[3] & min.cons.seg[,3] >= a[2])[1]
  geneMinConsOutRow   <- min.cons.seg[geneMinConsOutIndex,]
  
  # next, check how many extra genes we should look at
  startPosToConsiderIndex <- min(which(min.cons.seg[,1] %in% a[1] & min.cons.seg[,2]<= a[3] & min.cons.seg[,3] >= min.cons.seg[geneMinConsOutIndex,2]-basePairsEitherSide))
  endPosToConsiderIndex   <- max(which(min.cons.seg[,1] %in% a[1] & min.cons.seg[,2]<= min.cons.seg[geneMinConsOutIndex,3]+basePairsEitherSide & min.cons.seg[,3] >= a[2]))
  
  min.cons.seg          <- min.cons.seg[startPosToConsiderIndex:endPosToConsiderIndex,]
  min.cons.seg.gain     <- min.cons.seg.gain[startPosToConsiderIndex:endPosToConsiderIndex,]
  min.cons.seg.loss     <- min.cons.seg.loss[startPosToConsiderIndex:endPosToConsiderIndex,]
  perc.min.copy.gain    <- perc.min.copy.gain[startPosToConsiderIndex:endPosToConsiderIndex]
  perc.min.copy.loss    <- perc.min.copy.loss[startPosToConsiderIndex:endPosToConsiderIndex]
  min.cons.out[[1]]     <- min.cons.out[[1]][startPosToConsiderIndex:endPosToConsiderIndex,]
  min.cons.out[[2]]     <- min.cons.out[[2]][startPosToConsiderIndex:endPosToConsiderIndex,]
  min.cons.out[[3]]     <- min.cons.out[[3]][startPosToConsiderIndex:endPosToConsiderIndex,]
  min.cons.out[[4]]     <- min.cons.out[[4]][startPosToConsiderIndex:endPosToConsiderIndex,]
  min.cons.out[[5]]     <- min.cons.out[[5]][startPosToConsiderIndex:endPosToConsiderIndex,]
  min.cons.out[[6]]     <- min.cons.out[[6]][startPosToConsiderIndex:endPosToConsiderIndex,]
  
  #prepare data for plotting
  plot.pos <- min.cons.seg[,1:3]
  act.pos  <- min.cons.seg[,1:3]
  chr.size <- max(plot.pos[,3])
  chr.size.cum <- chr.size

  plot.gain.vect <- c(setNames(perc.min.copy.gain, plot.pos[,2]), setNames(perc.min.copy.gain, plot.pos[,3]))
  plot.gain.vect <- plot.gain.vect[order(as.numeric(names(plot.gain.vect)))]
  plot.loss.vect <- c(setNames(perc.min.copy.loss, plot.pos[,2]), setNames(perc.min.copy.loss, plot.pos[,3]))
  plot.loss.vect <- plot.loss.vect[order(as.numeric(names(plot.loss.vect)))]
  
  #if type == 'all' plot clonal and subclonal frequencies
  if(type == 'all'){
    freq.min.copy.loss_trunk     <-  apply(min.cons.out[[4]][,sample.names],1,sum,na.rm=T)
    perc.min.copy.loss_trunk     <-  freq.min.copy.loss_trunk/ncol(min.cons.out[[4]][,sample.names])
    freq.min.copy.gain_trunk     <-  apply(min.cons.out[[3]][,sample.names],1,sum,na.rm=T)
    perc.min.copy.gain_trunk     <-  freq.min.copy.gain_trunk/ncol(min.cons.out[[3]][,sample.names])
    
    plot.gain_trunk.vect <- c(setNames(perc.min.copy.gain_trunk, plot.pos[,2]), setNames(perc.min.copy.gain_trunk, plot.pos[,3]))
    plot.gain_trunk.vect <- plot.gain_trunk.vect[order(as.numeric(names(plot.gain_trunk.vect)))]
    plot.loss_trunk.vect <- c(setNames(perc.min.copy.loss_trunk, plot.pos[,2]), setNames(perc.min.copy.loss_trunk, plot.pos[,3]))
    plot.loss_trunk.vect <- plot.loss_trunk.vect[order(as.numeric(names(plot.loss_trunk.vect)))]
  }
  
  cols <- setNames(c('#fb6a4a','#bd0026','#9ecae1','#0571b0','#31a354','#006d2c','#e6550d','#a63603'), c('branch_gain','trunk_gain','branch_loss','trunk_loss','gistic_branch_gain','gistic_trunk_gain','gistic_branch_loss','gistic_trunk_loss'))
  
  geneIndex <- which(act.pos[,1] %in% a[1] & act.pos[,2]<= a[3] & act.pos[,3] >= a[2])[1]
  geneChrom <- a[1]
  
  geneIndexMinCons <- which(min.cons.seg.gain[,1] %in% a[1] & min.cons.seg.gain[,2]< a[3] & min.cons.seg.gain[,3] > a[2])[1]
  geneSpec         <- cbind(plot.loss.vect[as.character(plot.pos[geneIndex,2])], plot.gain.vect[as.character(plot.pos[geneIndex,2])])
  
  #open pdf
  pdfname <- paste0(pdf_name,'_',GeneOfInterest,'.pdf')
  pdf(pdfname, width=10, height=6, useDingbats = F)
  par(las=1)
  plot(perc.min.copy.gain, type='n', ylim=c(-1.1,1.1), xlim=c(min(act.pos[,2])/1e6,sum(chr.size)/1e6),main=main,ylab='Frequency aberrant', xlab=paste0('Chromosome ', geneChrom, ' (mb)'),axes=F)
  
  polygon(c(as.numeric(names(plot.gain.vect))/1e6, rev(as.numeric(names(plot.gain.vect))/1e6)), c(plot.gain.vect,rep(0,length(plot.gain.vect))), col=cols['branch_gain'],border=NA)
  polygon(c(as.numeric(names(plot.loss.vect))/1e6, rev(as.numeric(names(plot.loss.vect))/1e6)), c(plot.loss.vect,rep(0,length(plot.loss.vect))), col=cols['branch_loss'],border=NA)
  abline(h=0, col=1)
  #plot clonal and subclonal frequencies
  if(type == 'all'){
    polygon(c(as.numeric(names(plot.gain_trunk.vect))/1e6, rev(as.numeric(names(plot.gain_trunk.vect))/1e6)), c(plot.gain_trunk.vect,rep(0,length(plot.gain_trunk.vect))), col=cols['trunk_gain'],border=NA)
    polygon(c(as.numeric(names(plot.loss_trunk.vect))/1e6, rev(as.numeric(names(plot.loss_trunk.vect))/1e6)), c(plot.loss_trunk.vect,rep(0,length(plot.loss_trunk.vect))), col=cols['trunk_loss'],border=NA)
  }
  
  #annotate additional genes
  if(!is.null(add.gene.pos)) {
    for (x in add.gene.pos$symbol){
      gene.pos <- add.gene.pos[add.gene.pos[,1] == x,]
      ww       <- which(act.pos[,1] %in% as.numeric(sub('chr', '', gene.pos[2])) & act.pos[,2]<= as.numeric(gene.pos[4]) & act.pos[,3] >= as.numeric(gene.pos[3]))[1]
      y.pos    <- cbind(plot.loss.vect[as.character(plot.pos[ww,2])], plot.gain.vect[as.character(plot.pos[ww,2])])
      points(as.numeric(gene.pos[3])/1e6,y.pos[,2], pch = 18, col = '#660066')
    }
  }
  
  #add gene of interest with error
  points(as.numeric(rownames(geneSpec))/1e6,geneSpec[,2], pch = 16, col = '#00e6e6')
  arrowHeight <- as.numeric(geneSpec[,2])+0.05
  arrowPos <- as.numeric(rownames(geneSpec))/1e6
  arrows(x0=arrowPos,y0=1,x1=arrowPos,y1=arrowHeight,col='red',angle = 15,lwd=2)
  text(GeneOfInterest,x = arrowPos,y=1.05)
  
  #add line for quantiles derived from bootstrapping
  if(add_pval){
    abline(h=c(-sig.thresh.copy[1], sig.thresh.copy[2],0), lty=c(2,2,1), lwd=2, col=c(6,6,1))
  }

  #add line to seperate p and q arm
  cytobands     <- read.table(cytobands.file, header = F, stringsAsFactors = F)
  cytobands[,1] <- sub('chr', '', cytobands[,1])
  chr.cyto      <- cytobands[cytobands[,1] == geneChrom,]
  ww            <- min(grep('q', chr.cyto[,4]))
  split.pq      <- chr.cyto[ww, 2]
  if(split.pq > min(act.pos[,2]) & split.pq < max(act.pos[,3])){
    abline(v = split.pq/1e6, lty = 2)
    p.pos <- (split.pq - min(act.pos[,2])) / 2 + min(act.pos[,2])
    q.pos <- (max(act.pos[,3]) - split.pq) / 2 + split.pq
    text('p-arm', x = p.pos/1e6, y = -1)
    text('q-arm', x = q.pos/1e6, y = -1)
  }
  
  #add axis
  axis(2)
  axis(1, at=c(min(act.pos[,2])/1e6,sum(chr.size)/1e6), labels=c(min(act.pos[,2]),sum(chr.size)))
  box()
  dev.off()
  
  #combine data for output
  out.list <- list(plot.gain.vect, plot.gain_trunk.vect, plot.loss.vect, plot.loss_trunk.vect, sig.thresh.copy, sample.names)
  names(out.list) <- c('frequency.gain', 'frequency.gain_trunk', 'frequency.loss', 'frequency.loss_trunk', 'sig.thresh.copy', 'Samples')
  return(out.list)
}




##### function to get cn status of gene list in samples #####
findGene.cn <- function(genelist
                        , seg
                        , loss=log2(1.5/2)
                        , min.loss=log2(1.5/2)
                        , gain=log2(2.5/2)
                        , min.gain=log2(2.5/2)){
  
  #create empty sample_summary table
  sample_summary           <- matrix(0,0,19)
  colnames(sample_summary) <- c('SampleID','Gene','chr','start','stop','nregions','gained','lost','min_gained','min_lost','maxCN','minCN','gained.regions','lost.regions','allCN','allCNa','allCNb','allPloidy','allACF')
  
  #load gene positions
  genepos <- get.genepos(load.genepos=TRUE, return.genepos=TRUE,posfile = UCSCgenes.file)
  
  #get copy number of genes for each sample in seg
  samples <- unique(seg[,1])
  for(i in 1:length(samples)){
    sample.seg     <- seg[seg[,1] %in% samples[i],]
    
    #extract copy number for each gene in genelist
    for(j in 1:length(genelist)){
      pos <- tryCatch(get.genepos(genelist[j], genepos = genepos,posfile = UCSCgenes.file), error=function(err) return(NA))
      if(is.na(pos[1])){next}
      gene.seg   <- get.cn(sample.seg, pos[1], pos[2], pos[3], adjusted.seg=TRUE)
      gained     <- sum(log2(gene.seg[,6]/gene.seg[,9]) > gain)
      gained.regions <- gene.seg[log2(gene.seg[,6]/gene.seg[,9]) > gain,1]
      if(length(gained.regions) > 0) {
        gained.regions <- paste(gained.regions, collapse = '; ')
      } else {
        gained.regions <- NA
      }
      min.gained <- min.lost <- 0
      if(any(gained)){    
        min.gained <- sum(log2(gene.seg[,6]/gene.seg[,9]) > min.gain)
      }
      lost <- sum(log2(gene.seg[,6]/gene.seg[,9]) < loss)
      lost.regions <- gene.seg[log2(gene.seg[,6]/gene.seg[,9]) < loss, 1]
      if(length(lost.regions) > 0) {
        lost.regions <- paste(lost.regions, collapse = '; ')
      } else {
        lost.regions <- NA
      }
      if(any(lost)){    
        min.lost <- sum(log2(gene.seg[,6]/gene.seg[,9]) < min.loss)
      }
      if(any(c(gained,lost))){
        tmp <- cbind(samples[i], genelist[j],pos[1],pos[2],pos[3], length(sample_columns), gained, lost,  
                     min.gained, min.lost, max(gene.seg[,6]), min(gene.seg[,6]), gained.regions, lost.regions,
                     paste(paste(unique(gene.seg[,1]),gene.seg[,6],sep=':', collapse='; ')),
                     paste(paste(unique(gene.seg[,1]),gene.seg[,7],sep=':', collapse='; ')),
                     paste(paste(unique(gene.seg[,1]),gene.seg[,8],sep=':', collapse='; ')),
                     paste(paste(unique(gene.seg[,1]),signif(gene.seg[,9],2),sep=':', collapse='; ')),
                     paste(paste(unique(gene.seg[,1]),gene.seg[,10],sep=':', collapse='; ')))
        colnames(tmp)  <- c('SampleID','Gene','chr','start','stop','nregions','gained','lost','min_gained','min_lost','maxCN','minCN','gained.regions','lost.regions','allCN','allCNa','allCNb','allPloidy','allACF')
        sample_summary <- rbind(sample_summary,tmp)
      } 
    }
  }
  if(nrow(sample_summary) > 0){
    rownames(sample_summary) <- 1:nrow(sample_summary)
    sample_summary           <- as.data.frame(sample_summary, as.is=T)
    for(i in 3:10){
      sample_summary[,i] <- as.numeric(as.character(sample_summary[,i]))
    }
  }
  return(sample_summary)
}


##### function to test expression of a gene in samples where its gained/lost vs not gained/lost ###
test_cn.genes.expression <- function(freq.cn #(= output from freq_plot_gene.fn())
                                     , seg
                                     , expr.table
                                     , gain.threshold=log2(2.5/2)
                                     , loss.threshold=log2(1.5/2)
                                     , cancer.type
                                     , type = 'gain'
                                     , chr = 'chr1'
                                     , GeneOfInterest = 'MCL1'
                                     , plot.prefix) {
  #get genes that fall into amplified region and classify them as gained and lost
  if(type == 'gain') {
    sig.thresh.gain <- freq.cn[[5]][2]
    freq.gain       <- freq.cn[[1]]
    region.start <- min(as.numeric(names(freq.gain[freq.gain >= sig.thresh.gain])))
    region.stop  <- max(as.numeric(names(freq.gain[freq.gain >= sig.thresh.gain])))
  }
  if(type == 'loss') {
    sig.thresh.loss <- freq.cn[[5]][1]
    freq.loss       <- freq.cn[[3]]
    region.start <- min(as.numeric(names(freq.loss[freq.loss <= -1*sig.thresh.loss])))
    region.stop  <- max(as.numeric(names(freq.loss[freq.loss <= -1*sig.thresh.loss])))
  }
  
  region       <- GRanges(seqnames = chr, ranges = IRanges(start = region.start, end = region.stop))
  region.genes <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), region)
  region.genes <- merge(as.data.frame(region.genes), as.data.frame(org.Hs.egSYMBOL), by = 'gene_id')
  
  genelist <- unique(region.genes$symbol)
  cn.genes <- findGene.cn(genelist, seg, gain=gain.threshold, loss=loss.threshold)
  
  #compare gene expression between samples that have a gain/loss vs no gain/loss
  all.samples <- unique(seg$sample)
  
  output.table <- c()
  for (gene in as.character(unique(cn.genes$Gene))){
    #print(gene)
    expr.values <- expr.table[rownames(expr.table) %in% gene,]
    if(length(expr.values) == 0 | sum(expr.values) == 0) {
      next
    }
    expr.values <- expr.values[!is.na(expr.values)]
    sub         <- cn.genes[cn.genes$Gene == gene, ]
    if(type == 'gain'){
      region.column = 'gained.regions'
    }
    if(type == 'loss'){
      region.column = 'lost.regions'
    }
    samples.withEvent <- sapply(1:nrow(sub), function(x){
      sample  <- sub$SampleID[x]
      regions <- as.character(sub[x, region.column])
      regions <- unlist(strsplit(regions, '; '))
      return(paste(sample, regions, sep = ':'))
    })
    samples.withEvent    <- unlist(samples.withEvent)
    samples.withoutEvent <- all.samples[!all.samples%in%samples.withEvent]
    
    expr.withEvent    <- expr.values[names(expr.values)%in%samples.withEvent]
    expr.withoutEvent <- expr.values[names(expr.values)%in%samples.withoutEvent]
    if(sum(expr.withEvent) == 0 | sum(expr.withoutEvent) == 0){
      next
    }
    
    #wilcox-test
    p.value <- wilcox.test(expr.withEvent, expr.withoutEvent, alternative = 'greater')$p.value
    
    out <- rbind(data.frame(gene = gene, expression = expr.withEvent, group = type, pvalue = p.value, stringsAsFactors = F),
                 data.frame(gene = gene, expression = expr.withoutEvent, group = paste0('no', type), pvalue = p.value, stringsAsFactors = F))
    rownames(out) <- NULL
    output.table  <- rbind(output.table, out)
  }
  
  if(!is.null(output.table)) {
    #plot boxplot of significant genes
    plot.table <- output.table[output.table$pvalue < 0.05,]
    pl <- lapply(unique(plot.table$gene), function(x){
      data <- plot.table[plot.table$gene == x,]
      p <- ggplot(data, aes(x = group, y = expression, fill = group)) +
        geom_boxplot(outlier.size=-1) +
        geom_jitter(shape=16, position=position_jitter(0.2)) +
        stat_compare_means(method = 'wilcox.test', method.args = list(alternative = "less")) +
        theme_bw() +
        theme(legend.position = 'none') +
        ggtitle(x)
      return(p)
    })
    
    ml <- marrangeGrob(pl, nrow=2, ncol=2)
    pdf(paste0(plot.prefix, '.sigExpression.genes.', type, '.pdf'), width = 5, height = 6)
    ml
    dev.off()
    #ggsave(filename = paste0(plot.prefix, '.sigExpression.genes.', type, '.pdf'), plot = ml, width = 5, height = 6)
    
    #plot boxplot of gene of interest
    data <- output.table[output.table$gene == GeneOfInterest,]
    if(nrow(data) == 0) {
      print(paste0('Expression of ', GeneOfInterest, ' is 0 in all samples!'))
    } else {
      p <- ggplot(data, aes(x = group, y = expression, fill = group)) +
        geom_boxplot(outlier.size=-1) +
        geom_jitter(shape=16, position=position_jitter(0.2)) +
        stat_compare_means(method = 'wilcox.test', method.args = list(alternative = "less")) +
        theme_bw() +
        theme(legend.position = 'none') +
        ylab(paste0(GeneOfInterest, ' mRNA level (TPM)')) +
        ggtitle(GeneOfInterest)
      pdf(paste0(plot.prefix, '.expression.', type, '_', GeneOfInterest, '.pdf'), width = 4, height = 5, useDingbats = F)
      #tiff(file = paste0(plot.prefix, '.expression.', type, '_', GeneOfInterest, '.tiff'), width = 3200, height = 3600, units = "px", res = 800)
      plot(p)
      dev.off()
    }
    
    #plot barplot of -log(p.value)
    plot.pvalues            <- plot.table[!duplicated(plot.table$gene),]
    plot.pvalues$log.pvalue <- -log(plot.pvalues$pvalue)
    plot.pvalues            <- plot.pvalues[order(plot.pvalues$log.pvalue, decreasing = T),]
    plot.pvalues$colour     <- '#a6cee3'
    plot.pvalues$colour[plot.pvalues$gene == GeneOfInterest] <- '#e31a1c'
    pdf(paste0(plot.prefix, '.pvalues.sigExprGenes.pdf'), width = round(nrow(plot.pvalues)*0.15), height = 5)
    x <- barplot(plot.pvalues$log.pvalue, col = plot.pvalues$colour, 
                 main = paste0('Significant Genes ', cancer.type), ylab = '-log(pvalue)', xlab = 'genes')
    col <- sub('#a6cee3', 'black', plot.pvalues$colour)
    text(cex=0.7, x=x, y=-0.1, plot.pvalues$gene, xpd=TRUE, srt=45, adj = 1, col = col)
    dev.off()
    
    #pie chart with numbers of sig genes and not sig genes
    n.sigGenes    <- nrow(plot.pvalues)
    n.notsigGenes <- length(genelist) - n.sigGenes
    pdf(paste0(plot.prefix, '.pie.sigExprGenes.pdf'), width = 5, height = 5)
    pie(c(n.sigGenes, n.notsigGenes), 
        labels = c(paste0('significant (', n.sigGenes, ')'), paste0('not significant (', n.notsigGenes, ')')),
        main = cancer.type, 
        col = c('#009999', '#8f246b'))
    dev.off()
    
    #plot significant genes into frequency plot
    sig.genes    <- unique(plot.table$gene)
    sig.gene.pos <- region.genes[region.genes$symbol %in%sig.genes, c(7,2,3,4)]
    freq.cn      <- freq_plot_gene.fn(seg, GeneOfInterest = GeneOfInterest,
                                      main = cancer.type, 
                                      gain.threshold = gain.threshold,
                                      loss.threshold = loss.threshold,
                                      pdf_name = paste0(plot.prefix, '.cnFrequency.sigExprGenes'),
                                      add.gene.pos = sig.gene.pos)
    
    
    return(list(output.table, cn.genes))
  } else {
    print('No significany genes were found!')
    return(cn.genes)
  }
}



##### Heatmap plot with segment lengths, TP53 status and clonality #####

#get gene position
genepos       <- get.genepos(load.genepos=TRUE, return.genepos=TRUE,posfile = UCSCgenes.file)
gene.position <- get.genepos(gene,genepos)

#get cytoband
cytobands     <- read.table(cytoband.file, header = F, stringsAsFactors = F)
cytobands[,1] <- sub('chr', '', cytobands[,1])
ww            <- which(cytobands[,1] == gene.position[1] & cytobands[,2] <= gene.position[2] & cytobands[,3] >= gene.position[3])
cytoband      <- cytobands[ww,]
arm           <- substr(cytoband[4], start = 1, stop = 1)
tmp           <- cytobands[cytobands[,1] == gene.position[1],]
chr.arm.pos   <- data.frame(chr = unique(tmp[,1]), start = tmp[min(grep(arm, tmp[,4])),2], stop = tmp[max(grep(arm, tmp[,4])),3], arm = arm)

#get segments per patient
plot.table <- c()
for(sample in unique(seg$sample)) {
  sub.seg <- seg[seg$sample == sample,]
  ww      <- which(sub.seg$chr == chr.arm.pos$chr & sub.seg$startpos >= chr.arm.pos$start & sub.seg$endpos <= chr.arm.pos$stop)
  data    <- sub.seg[ww,]
  #data$cnTotal <- (data$nAraw + data$nBraw)
  data$cnTotal <- log2(data$cnTotal / data$Ploidy)
  data$cnTotal[data$cnTotal %in% '-Inf'] <- -3
  data$group   <- 'neutral'
  data$group[data$cnTotal > gain.threshold] <- 'gain'
  data$group[data$cnTotal < loss.threshold] <- 'loss'
  data$geneStatus <- 'FALSE'
  data$geneStatus[data$startpos <= gene.position[2] & data$endpos >= gene.position[3]] <- 'TRUE'
  
  plot.table <- rbind(plot.table, data)
}

#order by size of segment with gene of Interest
order.samples      <- plot.table[as.logical(plot.table$geneStatus),]
order.samples      <- order.samples[order(as.numeric(order.samples$endpos) - as.numeric(order.samples$startpos), decreasing = F),]
order.samples      <- order.samples[order(order.samples$group, decreasing = T),]
order.sample.names <- order.samples$sample
index.names        <- setNames(seq(0, length(order.sample.names)-1), order.sample.names)

#add clonality of gains
freqGains       <- findGene.cn(gene,seg,loss = loss.threshold,min.loss = loss.threshold, gain = gain.threshold, min.gain = gain.threshold)
freqGains       <- freqGains[freqGains$gained>=1,,drop=FALSE]
GainsHomo       <- freqGains[freqGains$nregions==freqGains$gained,,drop=FALSE]
GainsHetero     <- freqGains[freqGains$nregions!=freqGains$gained,,drop=FALSE]
clonalGains     <- as.character(GainsHomo$SampleID)
subclonalGains  <- as.character(GainsHetero$SampleID)
order.samples$clonality <- 'noGain'
order.samples$clonality[substr(order.samples$sample, 1, 8) %in% clonalGains]    <- 'clonalGain'
order.samples$clonality[substr(order.samples$sample, 1, 8) %in% subclonalGains] <- 'subclonalGain'
order.samples$clonality[order.samples$clonality == 'subclonalGain' & order.samples$group %in% c('loss', 'neutral')] <- 'noGain'

#add TP53 status
TP53mutants <- unique(luad.mutTableRegion$RegionID[luad.mutTableRegion$DriverMut == TRUE & luad.mutTableRegion$Hugo_Symbol == 'TP53'])
order.samples$TP53status <- 'wt'
order.samples$TP53status[order.samples$sample %in% TP53mutants] <- 'mut'
order.samples$TP53status <- factor(order.samples$TP53status, levels = c('wt', 'mut'))

#plot
plot.table$sample   <- factor(plot.table$sample, levels = order.sample.names)
plot.table$ystart   <- index.names[plot.table$sample]
plot.table$yend     <- index.names[plot.table$sample] + 1
plot.table$group    <- factor(plot.table$group, levels = c('gain', 'neutral', 'loss'))

plot.table$cn.group <- cut(plot.table$cnTotal, breaks = c(rev(seq(loss.threshold, floor(min(plot.table$cnTotal)), -0.5)), seq(gain.threshold, ceiling(max(plot.table$cnTotal)), 0.5)))
levels.cn.group     <- names(table(plot.table$cn.group)[table(plot.table$cn.group) != 0])
plot.table$cn.group <- factor(plot.table$cn.group, levels = levels.cn.group)
n.blue <- length(grep('-', levels(plot.table$cn.group))) - 1
n.red  <- length(levels(plot.table$cn.group)) - n.blue

order.samples$sample <- factor(order.samples$sample, levels = order.sample.names)
order.samples$size   <- as.numeric(order.samples$endpos) - as.numeric(order.samples$startpos)

heatmap <- ggplot() +
  geom_rect(data = plot.table, aes(xmin = startpos / 1000000, xmax = endpos / 1000000, ymin = ystart, ymax = yend, fill = cn.group), color = 'black', size = 0.1) + 
  scale_fill_manual(name = 'log(Total Copy Number)', values = c(rev(brewer.pal(n.blue + 2, "Blues")[-(1:2)]), 'lightgray', brewer.pal(n.red + 2, "Reds")[-(1:2)]), labels = c(paste0('loss ', levels.cn.group[1:n.blue]), 'neutral', paste0('gain ', levels.cn.group[(n.blue+2):length(levels.cn.group)]))) +
  geom_vline(xintercept = gene.position[2] / 1000000, linetype="dashed", color = "black", size=1) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab(paste0("Chromosome arm ", chr.arm.pos$chr, chr.arm.pos$arm, ' (mb)')) +
  ggtitle('Copy Number Satus of Segments around MCL1 in LUAD') +
  theme_bw() +
  theme(strip.text.y = element_blank()) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(panel.grid = element_blank()) 

bar <- ggplot(order.samples, aes(x = sample, y = size / 1000000)) + 
  geom_bar(stat = 'identity', fill = '#ffbf80') + 
  scale_y_reverse(limits = c(max(order.samples$size / 1000000) + 2, min(order.samples$size / 1000000) - 2), expand = c(0, 0)) +
  ylab(paste0(gene, '\nsegment size\n (mb)')) + xlab("") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'none') + 
  theme(panel.grid = element_blank()) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

tile_clonality <- ggplot(order.samples, aes(x = sample, y = 1, fill = clonality)) + 
  geom_tile() + 
  scale_fill_manual(name = '', values = c('#f1b6da', '#bababa', '#9970ab')) +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'none') + 
  theme(panel.grid = element_blank()) +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title = element_blank()) 

tile_TP53 <- ggplot(order.samples, aes(x = sample, y = 1, fill = TP53status)) + 
  geom_tile() + 
  scale_fill_manual(name = '', values = c('#b2df8a', '#33a02c')) +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = 'none') + 
  theme(panel.grid = element_blank()) +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title = element_blank()) 

#combine plots
heatmap        <- ggplotGrob(heatmap + theme(plot.margin = unit(c(1, 0.5, 0.5, 0), "cm")))
tile_clonality <- ggplotGrob(tile_clonality + theme(plot.margin = unit(c(1, 0.05, 0.5, 0.05), "cm")))
tile_TP53      <- ggplotGrob(tile_TP53 + theme(plot.margin = unit(c(1, 0.05, 0.5, 0.05), "cm")))
bar            <- ggplotGrob(bar + theme(plot.margin = unit(c(1, 0.1, 0.5, 0.5), "cm")))

heatmap$heights        <- unit.pmax(heatmap$heights, tile_clonality$heights, tile_TP53$heights, bar$heights)
tile_clonality$heights <- unit.pmax(heatmap$heights, tile_clonality$heights, tile_TP53$heights, bar$heights)
tile_TP53$heights      <- unit.pmax(heatmap$heights, tile_clonality$heights, tile_TP53$heights, bar$heights)
bar$heights            <- unit.pmax(heatmap$heights, tile_clonality$heights, tile_TP53$heights, bar$heights)

combined.list        <- list(bar, tile_clonality, tile_TP53, heatmap)
combined.plot.layout <- matrix(c(rep(1, 3), rep(2,1), rep(3,1), rep(4,20)), nrow = 1)
combined.plot        <- arrangeGrob(grobs = combined.list, ncol = dim(combined.plot.layout)[2], nrow = 1, layout_matrix = combined.plot.layout)

grid.draw(combined.plot)


##### function to plot heatmap with copy number and mutation information of common drivers #####
plot_heatmap.cn.mutations <- function(genes
                                      , geneOfInterest = 'MCL1'
                                      , seg
                                      , mut.table
                                      , gain.threshold=log2(2.5/2)
                                      , loss.threshold=log2(1.5/2)
                                      , title
                                      , pdf.name) {
  
  #only use samples with copy number and mutation information
  all.regions  <- intersect(unique(seg$RegionID), mut.table$RegionID)
  
  cn.plot.table <- c()
  mut.plot.table <- c()
  for (gene in genes){
    #get gene position
    genepos       <- get.genepos(load.genepos=TRUE, return.genepos=TRUE,posfile = UCSCgenes.file)
    gene.position <- get.genepos(gene,genepos)
    
    #get copy number status
    ww      <- which(seg$chr == gene.position[1] & seg$startpos <= gene.position[2] & seg$endpos >= gene.position[3])
    data.cn <- seg[ww,]
    data.cn$cnTotal <- (data.cn$nAraw + data.cn$nBraw)
    data.cn$cnTotal <- log2(data.cn$cnTotal / data.cn$Ploidy)
    data.cn$group   <- 'neutral'
    data.cn$group[data.cn$cnTotal > gain.threshold] <- 'gain'
    data.cn$group[data.cn$cnTotal < loss.threshold] <- 'loss'
    data.cn       <- data.frame(data.cn[,c('RegionID', 'group')], gene = gene, stringsAsFactors = F)
    cn.plot.table <- rbind(cn.plot.table, data.cn)
    
    #get mutation status
    ww       <- grep(paste0('^', gene, '$'), mut.table$Hugo_Symbol)
    data.mut <- mut.table[ww, c('RegionID', 'regionMutID', 'mutation_id', 'Hugo_Symbol', 'combTiming', 'exonic.func')]
    data.mut <- data.mut[!duplicated(data.mut$RegionID),]
    mut.plot.table <- rbind(mut.plot.table, data.mut)
  }
  
  #use samples that are present in seg and mut.table
  mut.plot.table           <- mut.plot.table[mut.plot.table$RegionID %in% all.regions,] 
  if(length(unique(mut.plot.table$RegionID)) != length(all.regions)) {
    add.samples <- all.regions[!all.regions %in% unique(mut.plot.table$RegionID)]
    add.df      <- data.frame(RegionID = add.samples, regionMutID = NA, mutation_id = NA, Hugo_Symbol = NA,
                              combTiming = NA, exonic.func = NA, stringsAsFactors = F)
    mut.plot.table <- rbind(mut.plot.table, add.df)
  }
  colnames(mut.plot.table) <- sub('Hugo_Symbol', 'gene', colnames(mut.plot.table))
  
  cn.plot.table <- cn.plot.table[cn.plot.table$RegionID %in% all.regions,]
  if(length(unique(cn.plot.table$RegionID)) != length(all.regions)) {
    add.samples <- all.regions[!all.regions %in% unique(cn.plot.table$RegionID)]
    add.df      <- data.frame(RegionID = add.samples, group = NA, gene = NA, stringsAsFactors = F)
    cn.plot.table <- rbind(cn.plot.table, add.df)
  }
  
  if(!is.null(geneOfInterest)) {
    sub                              <- cn.plot.table[cn.plot.table$gene == geneOfInterest,]
    mut.plot.table$geneOfInterest.cn <- sub$group[match(mut.plot.table$RegionID, sub$RegionID)]
    mut.plot.table <- mut.plot.table[order(mut.plot.table$geneOfInterest.cn),]
    order.samples  <- unique(mut.plot.table$RegionID)
  } else {
    order.samples <- all.regions
  }
  
  cn.plot.table$RegionID <- factor(cn.plot.table$RegionID, levels = order.samples)
  cn.plot.table$gene     <- factor(cn.plot.table$gene, levels = genes)
  cn.plot.table$group    <- factor(cn.plot.table$group, levels = c('gain', 'neutral', 'loss'))
  
  mut.plot.table$RegionID    <- factor(mut.plot.table$RegionID, levels = order.samples)
  mut.plot.table$gene        <- factor(mut.plot.table$gene, levels = genes)
  mut.plot.table$combTiming  <- factor(mut.plot.table$combTiming, levels = c('early', 'late', 'subclonal', 'Unknown'))
  mut.plot.table$exonic.func <- as.factor(mut.plot.table$exonic.func)
  #mut.plot.table$exonic.func[is.na(mut.plot.table$exonic.func)] <- 'unknown'
  
  #plot with timing of mutations
  p <- ggplot() +
    geom_tile(data = cn.plot.table, aes(x=RegionID, y = gene, fill = group), colour = 'black') +
    scale_fill_manual(name = 'Copy Number Status', values = alpha(c('gain' = '#fb8072', 'neutral' = '#d9e6f2', 'loss' = '#80b1d3'), 0.5), labels = c('gain', 'neutral', 'loss')) +
    geom_point(data = mut.plot.table[!is.na(mut.plot.table$gene),], aes(x = RegionID, y = gene, colour = exonic.func, shape = combTiming)) +
    #scale_color_manual(name = 'Timing', values = c('early' = '#d147a3', 'late' = '#732673', 'subclonal' = '#009999', 'Unknown' = '#666633'), labels = c('early', 'late', 'subclonal', 'unknown')) +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank()) +
    xlab('Samples') + ylab('Genes') + ggtitle(title)
  
  pdf(pdf.name, width = 15, height = 5, useDingbats = F)
  plot(p)
  dev.off()
  
  return(list(mut.plot.table, cn.plot.table))
}





