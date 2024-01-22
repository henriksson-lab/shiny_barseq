library(umap)
library(stringr)
library(plotly)
library(reshape2) #was reshape
library(sqldf)

### To generate a count matrix:
### /corgi/otherdataset/ellenbushell/barseq_pools$ python ../count_barseq.py slowhires_2023dec


dir_fancyplot <- "/corgi/otherdataset/ellenbushell/plots_barseq_fancy"
genedesc <- read.csv("/corgi/otherdataset/ellenbushell/gene_description.csv",sep="\t")[,c(1,2)]



list_pools <- list.files("/corgi/otherdataset/ellenbushell/barseq_pools")
list_pools <- c(
  
  "slowhires_2023dec",

  "EB_minipool2",
  
  "EB_barseq_slowpool_1","EB_barseq_slowpool_2",
  "merged_barseq_slowpool_1_and_2",
  
  "EB_priming_barseqpool1","EB_priming_barseqpool3","EB_priming_barseqpool4",
  "merged_EB_priming_barseqpool",
  
  "EB_priming_Candidatepool1","EB_priming_Candidatepool2",
  "merged_EB_priming_Candidatepool"
)

#list_pools <- "slowhires_2023dec"

################################################################################
###################### merging of pools ########################################
################################################################################

##################### merging of slowpool 1+2
if(FALSE){
  ### Merge

  dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_barseq_slowpool_1/counts.csv"))
  dat <- dat[!duplicated(dat$X),] #apparently there are dups!
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat1 <- dat
  
  dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_barseq_slowpool_2/counts.csv"))
  dat <- dat[!duplicated(dat$X),] #apparently there are dups!
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat2 <- dat

  dat12 <- cbind(dat1,dat2[rownames(dat1),])
  write.csv(dat12, "/corgi/otherdataset/ellenbushell/barseq_pools/merged_barseq_slowpool_1_and_2/counts.csv")
  #dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/merged_barseq_slowpool_1_and_2/counts.csv"))
  #head(dat)  
  
  samplemeta1 <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_barseq_slowpool_1/sampleinfo.txt"),sep="\t")
  samplemeta2 <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_barseq_slowpool_2/sampleinfo.txt"),sep="\t")

  samplemeta1$User.ID <- str_replace_all(samplemeta1$User.ID,"_m1","_m11")
  samplemeta1$User.ID <- str_replace_all(samplemeta1$User.ID,"_m2","_m12")
  samplemeta1$User.ID <- str_replace_all(samplemeta1$User.ID,"_m3","_m13")
  
  samplemeta2$User.ID <- str_replace_all(samplemeta2$User.ID,"_m2","_m22")
  samplemeta2$User.ID <- str_replace_all(samplemeta2$User.ID,"_m1","_m21")
  samplemeta2$User.ID <- str_replace_all(samplemeta2$User.ID,"_m3","_m23")

  samplemeta1$User.ID <- str_replace_all(samplemeta1$User.ID,"poolSlow1","poolSlow")
  samplemeta2$User.ID <- str_replace_all(samplemeta2$User.ID,"poolSlow2","poolSlow")
  
  head(samplemeta2)
  
  write.table(
    rbind(samplemeta1,samplemeta2),
    "/corgi/otherdataset/ellenbushell/barseq_pools/merged_barseq_slowpool_1_and_2/sampleinfo.txt",
    row.names = FALSE,sep="\t"
  )
}


##################### merging of primingpool
if(FALSE){
  ### Merge
  
  dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_Candidatepool1/counts.csv"))
  dat <- dat[!duplicated(dat$X),] #apparently there are dups!
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat1 <- dat
  
  dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_Candidatepool2/counts.csv"))
  dat <- dat[!duplicated(dat$X),] #apparently there are dups!
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat2 <- dat
  
  merged_dat <- cbind(dat1,dat2[rownames(dat1),])
  write.csv(merged_dat, "/corgi/otherdataset/ellenbushell/barseq_pools/merged_EB_priming_Candidatepool/counts.csv")

  samplemeta1 <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_Candidatepool1/sampleinfo.txt"),sep="\t")
  samplemeta2 <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_Candidatepool2/sampleinfo.txt"),sep="\t")
  
  samplemeta1$User.ID <- str_replace_all(samplemeta1$User.ID,"_m1","_m11")
  samplemeta1$User.ID <- str_replace_all(samplemeta1$User.ID,"_m2","_m21")
  samplemeta1$User.ID <- str_replace_all(samplemeta1$User.ID,"_m3","_m31")
  
  samplemeta2$User.ID <- str_replace_all(samplemeta2$User.ID,"_m1","_m12")
  samplemeta2$User.ID <- str_replace_all(samplemeta2$User.ID,"_m2","_m22")
  samplemeta2$User.ID <- str_replace_all(samplemeta2$User.ID,"_m3","_m32")
  
  write.table(
    rbind(samplemeta1,samplemeta2),
    "/corgi/otherdataset/ellenbushell/barseq_pools/merged_EB_priming_Candidatepool/sampleinfo.txt",
    row.names = FALSE,sep="\t"
  )
}

##################### merging of EB_priming_barseqpool1+3+4
if(FALSE){
  ### Merge
  
  dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_barseqpool1/counts.csv"))
  dat <- dat[!duplicated(dat$X),] #apparently there are dups!
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat1 <- dat
  
  dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_barseqpool3/counts.csv"))
  dat <- dat[!duplicated(dat$X),] #apparently there are dups!
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat3 <- dat

  dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_barseqpool4/counts.csv"))
  dat <- dat[!duplicated(dat$X),] #apparently there are dups!
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat4 <- dat
  
  merged_dat <- cbind(dat1,dat2[rownames(dat1),],dat3[rownames(dat1),])
  write.csv(merged_dat, "/corgi/otherdataset/ellenbushell/barseq_pools/merged_EB_priming_barseqpool/counts.csv")
  #dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/merged_barseq_slowpool_1_and_2/counts.csv"))
  #head(dat)  
  
  samplemeta1 <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_barseqpool1/sampleinfo.txt"),sep="\t")
  samplemeta3 <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_barseqpool3/sampleinfo.txt"),sep="\t")
  samplemeta4 <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_barseqpool4/sampleinfo.txt"),sep="\t")
  
  samplemeta1$User.ID <- str_replace_all(samplemeta1$User.ID,"_m1","_m11")
  samplemeta1$User.ID <- str_replace_all(samplemeta1$User.ID,"_m2","_m21")
  samplemeta1$User.ID <- str_replace_all(samplemeta1$User.ID,"_m3","_m31")

  samplemeta3$User.ID <- str_replace_all(samplemeta3$User.ID,"_m1","_m13")
  samplemeta3$User.ID <- str_replace_all(samplemeta3$User.ID,"_m2","_m23")
  samplemeta3$User.ID <- str_replace_all(samplemeta3$User.ID,"_m3","_m33")

  samplemeta4$User.ID <- str_replace_all(samplemeta4$User.ID,"_m1","_m14")
  samplemeta4$User.ID <- str_replace_all(samplemeta4$User.ID,"_m2","_m24")
  samplemeta4$User.ID <- str_replace_all(samplemeta4$User.ID,"_m3","_m34")
  
  #samplemeta1$User.ID <- str_replace_all(samplemeta1$User.ID,"poolSlow1","poolSlow")
  #samplemeta2$User.ID <- str_replace_all(samplemeta2$User.ID,"poolSlow2","poolSlow")
  
  #head(samplemeta2)
  
  write.table(
    rbind(samplemeta1, samplemeta3, samplemeta4),
    "/corgi/otherdataset/ellenbushell/barseq_pools/merged_EB_priming_barseqpool/sampleinfo.txt",
    row.names = FALSE,sep="\t"
  )
}






################################################################################
###################### all QC ##################################################
################################################################################



####### Sequencing depth per sample
for(current_pool in list_pools){
  print("Sequencing depth per sample")
  print(current_pool)

  dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/",current_pool,"/counts.csv"))
  dat <- dat[!duplicated(dat$X),] #apparently there are dups!
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  print(mean(as.matrix(dat)))
}  

###### Draw all basic RGR curves
all_samplemeta <- list()
for(current_pool in list_pools){
  print("Draw all basic RGR curves")
  print(current_pool)
  
  control_gene_file <- file.path("/corgi/otherdataset/ellenbushell/barseq_pools",current_pool,"/controls.csv")
  if(file.exists(control_gene_file)){
    print("using control gene list")
    control_genes <- read.csv(control_gene_file)
    control_genes$reason <- "control"
  } else {
    print("using default control genes")
    control_genes <- read.csv("/corgi/otherdataset/ellenbushell/dispensible_controls_for_barseq.csv",sep="\t")
  }
  
  
  #####################################
  
  #Only retain main mutants (some have multiple barcodes for them)
  dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/",current_pool,"/counts.csv"))
  dat <- dat[!duplicated(dat$X),] #apparently there are dups!
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]

  
  used_bc_file <- file.path("/corgi/otherdataset/ellenbushell/barseq_pools",current_pool,"/used_bc.csv")  #needs renaming to be used
  
  if(file.exists(used_bc_file)){
    print("using bc list")
    used_bc <- unique(c(control_genes$geneid, read.csv(used_bc_file)[,1]))
  } else {
    print("not using any bc list")
    #stop("this is no longer allowed")
    used_bc <- c()
  }
  
  ##################### How many reads on the expected BCs?
  
  bc_stats <- as.data.frame(rowSums(dat[-c(1:2),]))
  colnames(bc_stats) <- "sum"
  bc_stats$used <- FALSE
  if(length(used_bc)>0){
    bc_stats[used_bc,]$used <- TRUE
  }
  bc_stats <- bc_stats[order(bc_stats$sum, decreasing = TRUE),]
  head(bc_stats)
  bc_stats$log_sum <- log10(1+bc_stats$sum)
  
  bc_stats$x <- 1:nrow(bc_stats)
  ggplot(bc_stats, aes(x, log_sum, color=used)) + geom_point() +
    xlab("sgRNA index, sorted by abundance") + ylab("log10 abundance")
  ggsave(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "BC abundance.pdf"))) 
  write.csv(bc_stats, file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "BC abundance.csv"))) 
  
  
  ##################### Read sample meta
  samplemeta <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools",current_pool,"sampleinfo.txt"),sep="\t")
  colnames(samplemeta)[1:2] <- c("sampleid","barseqid")#,"k_reads","q30")

  samplemeta$primed <- "up"
  samplemeta$genotype <- "ug"

  samplemeta$genotype[grep("_RAG1KO_",samplemeta$barseqid)] <- "RAG1KO"
  samplemeta$genotype[grep("_Rag1_",samplemeta$barseqid)] <- "RAG1KO"
  samplemeta$genotype[grep("_BL6_",samplemeta$barseqid)] <- "BL6"
  
  samplemeta$genotype[grep("_IFNAR",samplemeta$barseqid)] <- "IFNAR"
  samplemeta$genotype[grep("_WT",samplemeta$barseqid)] <- "WT"
  
  samplemeta$primed[grep("_NP_",samplemeta$barseqid)] <- "NP"
  samplemeta$primed[grep("_P_",samplemeta$barseqid)] <- "P"
  
  samplemeta$pool <- str_split_fixed(str_split_fixed(samplemeta$barseqid,"pool",2)[,2],"_",2)[,1]
  samplemeta$pool[samplemeta$pool==""] <- "uP"
  
  samplemeta$day <- str_split_fixed(str_split_fixed(samplemeta$barseqid,"_d",2)[,2],"_",2)[,1]
  samplemeta$mouse <- str_split_fixed(str_split_fixed(samplemeta$barseqid,"_m",2)[,2],"_",2)[,1]
  
  samplemeta$is_input <- FALSE
  samplemeta$is_input[grep("_input",samplemeta$barseqid)] <- TRUE
  
  
  #Sum up reads for R1 and R2, for each sample
  cnt <- matrix(nrow = nrow(dat), ncol=nrow(samplemeta))
  rownames(cnt) <- rownames(dat)
  colnames(cnt) <- samplemeta$sampleid
  for(i in 1:nrow(samplemeta)){
    cnt[,i] <- rowSums(dat[,str_starts(colnames(dat), sprintf("%s.R",samplemeta$sampleid[i])), drop=FALSE])
  }

  ## Exclude bad samples
  exclude_sample_file <- file.path("/corgi/otherdataset/ellenbushell/barseq_pools/",current_pool,"/exclude_samples.csv")
  if(file.exists(exclude_sample_file)){
    exclude_samples <- read.csv(exclude_sample_file)[,1]
    keep <- !(samplemeta$sampleid %in% exclude_samples)
    samplemeta <- samplemeta[keep,]
    cnt <- cnt[,keep]
    print(paste("excluding samples"))
    print(exclude_samples)
  }
  
  
  samplemeta$counts_other <- cnt["_other",]
  samplemeta$counts_err <- cnt["_err",]
  samplemeta$counts_tot <- colSums(cnt)
  
  
  ##################### Only retain expected BCs, if list of those present
  if(length(used_bc)>0){
    print("only keeping BCs listed")
    cnt <- cnt[rownames(cnt) %in% used_bc,]
  } else {
    cnt <- cnt[-c(1:2),]  #remove first two at least
    
    #remove some more? keep top N barcodes?
    cnt <- cnt[rowSums(cnt)>0,]
  }
  
  #Compute error statistics
  samplemeta$perc_other <- samplemeta$counts_other / samplemeta$counts_tot
  samplemeta$perc_err <- samplemeta$counts_err / samplemeta$counts_tot
  samplemeta$perc_ok <- 1-samplemeta$perc_other - samplemeta$perc_err
  samplemeta$counts_ok <- colSums(cnt)
  
  pdf(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "perlib_counts_ok.pdf"))) 
  plot(samplemeta$counts_ok)
  dev.off()
  
  pdf(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "perlib_perc_ok.pdf"))) 
  plot(samplemeta$perc_ok)
  dev.off()
  
  #Remove weird sample
  keep <- samplemeta$barseqid!="Reads_for_unknown_index"
  cnt <- cnt[,keep]
  samplemeta <- samplemeta[keep,]
  
  
  #Remove samples with poor coverage				6666666666 nasty cutoff! should be %
  samplemeta$sample_ok <- samplemeta$counts_ok > 3000   
  write.csv(samplemeta,file.path("/corgi/otherdataset/ellenbushell/filtered_barseq_counts",paste(sep="",current_pool, ".all.samplemeta.csv")))
  write.csv(cnt,file.path("/corgi/otherdataset/ellenbushell/filtered_barseq_counts",paste(sep="",current_pool, ".all.counts.csv")))
  cnt <- cnt[,samplemeta$sample_ok]
  samplemeta <- samplemeta[samplemeta$sample_ok,]
  sum(!samplemeta$sample_ok)
  mean(samplemeta$sample_ok)
  
  #Store reduced count tables
  write.csv(samplemeta,file.path("/corgi/otherdataset/ellenbushell/filtered_barseq_counts",paste(sep="",current_pool, ".filtered.samplemeta.csv")))
  write.csv(cnt,file.path("/corgi/otherdataset/ellenbushell/filtered_barseq_counts",paste(sep="",current_pool, ".filtered.counts.csv")))
  #write.csv(cnt,file.path("/corgi/otherdataset/ellenbushell/filtered_barseq_counts",paste(current_pool, "counts.csv")))
  
  
  
  
  ######################### Normalization #####################
  
  #Normalize each sample
  cnt <- cnt[!str_starts(rownames(cnt),"_"),]
  for(i in 1:ncol(cnt)){
    cnt[,i] <- cnt[,i]/sum(cnt[,i])
  }
  
  #Heatmap, gene vs sample
  heatmap_cnt <- cnt
  colnames(heatmap_cnt) <- samplemeta$barseqid
  ht <- ComplexHeatmap::Heatmap(heatmap_cnt,
                                column_names_gp = grid::gpar(fontsize = 6),
                                row_names_gp = grid::gpar(fontsize = 6))
  pdf(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "heatmap.pdf")))
  ComplexHeatmap::draw(ht)
  dev.off()
  
  ht <- ComplexHeatmap::Heatmap(cor(heatmap_cnt),
                                column_names_gp = grid::gpar(fontsize = 6),
                                row_names_gp = grid::gpar(fontsize = 6))
  pdf(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "corr.pdf")))
  ComplexHeatmap::draw(ht)
  dev.off()
  #ggplot(thecorr, aes(Var1,Var2,fill=value)) + geom_tile()
  #ggcorrplot::ggcorrplot(cor(heatmap_cnt),tl.cex = 5)
  #ggsave(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "corr.pdf")))
  
  #Heatmap, gene vs sample
  heatmap_cnt <- log10(cnt+1e-4)
  colnames(heatmap_cnt) <- samplemeta$barseqid
  ht <- ComplexHeatmap::Heatmap(heatmap_cnt,
                                column_names_gp = grid::gpar(fontsize = 6),
                                row_names_gp = grid::gpar(fontsize = 6))
  pdf(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "heatmapLOG.pdf")))
  ComplexHeatmap::draw(ht)
  dev.off()
  
  ggcorrplot::ggcorrplot(cor(heatmap_cnt),tl.cex = 5)
  ggsave(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "corrLOG.pdf")))
  
    
  #Give each mouse an ID
  samplemeta$mouse_ref <- sprintf("%s_%s_%s", samplemeta$mouse, samplemeta$genotype, samplemeta$primed)
  samplemeta$mouse_ref[samplemeta$is_input] <- ""
  

  
  ################ Note: last samples have no barcodes! because new?
  
  #Perform UMAP to compare
  umap.settings <- umap.defaults
  umap.settings$n_neighbors <- min(umap.settings$n_neighbors, ncol(cnt))
  cnt.umap <- umap(t(cnt), config=umap.settings)
  samplemeta$umap1 <- cnt.umap$layout[,1]
  samplemeta$umap2 <- cnt.umap$layout[,2]
  
  p1 <- ggplot(samplemeta, aes(umap1,umap2,color=mouse_ref))+geom_point()
  p2 <- ggplot(samplemeta, aes(umap1,umap2,color=day))+geom_point()
  p3 <- ggplot(samplemeta, aes(umap1,umap2,color=mouse))+geom_point()
  p4 <- ggplot(samplemeta, aes(umap1,umap2,color=is_input))+geom_point()
  p5 <- ggplot(samplemeta, aes(umap1,umap2,color=genotype))+geom_point()
  p6 <- ggplot(samplemeta, aes(umap1,umap2,color=primed))+geom_point()
  p7 <- ggplot(samplemeta, aes(umap1,umap2,color=pool))+geom_point()
  p1/p2|p3/p4|p5/p6|p7
  ggsave(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "umap.pdf"))) 
  
  
  ########## Interactive plot
  fig <- plot_ly(type = 'scatter', mode = 'markers') 
  fig <- fig %>%
    add_trace(
      x = samplemeta$umap1, 
      y = samplemeta$umap2, 
      text= samplemeta$barseqid,
      hoverinfo = 'text',
      marker = list(color='green'),
      showlegend = F
    )
  fig
  
  htmlwidgets::saveWidget(fig, file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "umap.html")), selfcontained = F, libdir = "lib")


  all_samplemeta[[current_pool]] <- samplemeta
}

saveRDS(all_samplemeta, "/corgi/websites/shinybarseq/samplemeta.rds")













################################################################################
###################### all growth curves #######################################
################################################################################


all_timecourses <- list()
for(current_pool in list_pools){
  print("All growth curves")
  print(current_pool)
  
  control_gene_file <- file.path("/corgi/otherdataset/ellenbushell/barseq_pools",current_pool,"/controls.csv")
  if(file.exists(control_gene_file)){
    print("using control gene list")
    control_genes <- read.csv(control_gene_file)
    control_genes$reason <- "control"
  } else {
    print("using default control genes")
    control_genes <- read.csv("/corgi/otherdataset/ellenbushell/dispensible_controls_for_barseq.csv",sep="\t")
  }
  
  samplemeta <- read.csv(file.path("/corgi/otherdataset/ellenbushell/filtered_barseq_counts",paste(sep="",current_pool, ".filtered.samplemeta.csv")), row.names = "X")
  cnt <- read.csv(file.path("/corgi/otherdataset/ellenbushell/filtered_barseq_counts",paste(sep="",current_pool, ".filtered.counts.csv")), row.names = "X")
  
  #cnt <- cnt[rowMeans(cnt)>100,]  #original
  cnt <- cnt[rowMeans(cnt)>100,]
  print(dim(cnt))
  
  #Ensure days in order
  samplemeta$day_int <- as.integer(samplemeta$day)  
  cnt <- cnt[,order(samplemeta$day_int)]
  samplemeta <- samplemeta[order(samplemeta$day_int),]
  
  #Only keep sensible samples; no inputs
  keep <- !samplemeta$is_input & !is.na(samplemeta$day) & samplemeta$genotype!="ug"
  samplemeta <- samplemeta[keep,]
  cnt <- cnt[,keep]
  print(dim(cnt))
  
  #Only process samples if anything remains  
  #if(ncol(cnt)*nrow(cnt)>0){
  
  #Normalize each sample based on total count --- pseudo count added
  cnt.reg <- cnt+1
  for(i in 1:ncol(cnt.reg)){
    cnt.reg[,i] <- cnt.reg[,i]/sum(cnt.reg[,i])
  }
  
  #Normalize each sample based on total count
  for(i in 1:ncol(cnt)){
    cnt[,i] <- cnt[,i]/sum(cnt[,i])
  }
  
  #Give each mouse an ID
  samplemeta$mouse_ref <- sprintf("%s_%s_%s_%s", samplemeta$pool, samplemeta$mouse, samplemeta$genotype, samplemeta$primed)
  samplemeta$mouse_ref[samplemeta$is_input] <- ""
  
  #Give each condition an ID    
  samplemeta$condition_ref <- sprintf("%s_%s_%s", samplemeta$pool, samplemeta$genotype, samplemeta$primed)
  
  
  
  #### Compute GR for each strain
  all_gr <- NULL
  alldays <- unique(samplemeta$day_int)
  for(cur_mr in unique(samplemeta$mouse_ref)){
    
    #For each day
    for(curday in alldays[-1]){
      
      #Figure out which samples to compare. Ignore cases of missing samples        
      this_samp <- which(samplemeta$mouse_ref==cur_mr & samplemeta$day_int==curday)
      prev_samp <- which(samplemeta$mouse_ref==cur_mr & samplemeta$day_int==curday-1)
      if(length(this_samp)==1 & length(prev_samp)==1){
        
        all_gr <- rbind(
          all_gr,
          data.frame(
            mouse_ref=cur_mr,
            gene=rownames(cnt),
            day_int=curday,
            gr=cnt.reg[,this_samp]/cnt.reg[,prev_samp]
          )
        )
      }
    }  
  }
  
  if(FALSE){
    all_gr$mouse_ref
    all_gr[all_gr$mouse_ref=="SLOW_3_BL6_P" & all_gr$gene=="PBANKA_100360",]
    all_gr[all_gr$gene=="PBANKA_100360",]
  }
  
  
  #### Compute average GR for the control genes. exclude control genes if they ever go out of bounds
  control_gr <- all_gr[all_gr$gene %in% control_genes$geneid,]
  invalid_control_genes <- unique(control_gr$gene[is.nan(control_gr$gr) | is.infinite(control_gr$gr)])
  control_gr <- control_gr[!(control_gr$gene %in% invalid_control_genes),]
  control_gr_average <- sqldf("select mouse_ref, day_int, avg(gr) as control_gr from control_gr group by mouse_ref, day_int")
  
  #### Compute RGR, per mouse
  all_gr_meta <- merge(all_gr,samplemeta[,c("barseqid","sampleid","day_int","mouse_ref","condition_ref")])
  all_gr_meta <- merge(all_gr_meta,control_gr_average)
  all_gr_meta$rgr <- all_gr_meta$gr / all_gr_meta$control_gr
  
  #### Compute RGR, per condition and time point
  all_gr_meta_percond <- sqldf("select day_int, gene, avg(rgr) as rgr, condition_ref from all_gr_meta group by day_int, gene, condition_ref")
  
  
  
  #### What is the variance?
  all_gr_meta_percond$log_rgr <- log2(all_gr_meta_percond$rgr)
  ggplot(all_gr_meta_percond, aes(day_int, log_rgr, group=day_int)) + geom_violin(fill="blue") + geom_jitter(size=0.1)
  ggsave(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "violin_rgr_percond.pdf")))
  
  all_gr_meta$log_rgr <- log2(all_gr_meta$rgr)
  ggplot(all_gr_meta, aes(day_int, log_rgr, group=day_int)) + geom_violin(fill="gray") + geom_jitter(size=0.1)
  ggsave(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "violin_rgr_permouse.pdf")))
  
  #### Do all statistical tests
  totest <- merge(all_gr_meta_percond, samplemeta[,c("condition_ref","day_int","genotype","primed")])
  totest$cond <- paste(sep="_",totest$primed, totest$genotype)
  pvalcalc <- NULL
  for(cur_cond in unique(totest$cond)){
    sub_totest <- totest[totest$cond==cur_cond,]
    for(cur_gene in unique(sub_totest$gene)){
      onetest <- t.test(sub_totest$log_rgr[sub_totest$gene==cur_gene])
      pvalcalc <- rbind(pvalcalc, data.frame(cond=cur_cond, gene=cur_gene, p=onetest$p.value, fc=onetest$estimate))
    }
  }
  
  
  
  
  # 
  # #### Compare all conditions
  # #"plot delta RGR bwtween Primed Rags and Primed BL6."
  # for(cond1 in unique(totest$cond)){
  #   for(cond2 in unique(totest$cond)){
  #     if(cond1!=cond2){
  #       compared_cond <- merge(
  #         data.frame(
  #           gene=pvalcalc[pvalcalc$cond==cond1,]$gene,
  #           p1=-log10(pvalcalc[pvalcalc$cond==cond1,]$p),
  #           fc1=pvalcalc[pvalcalc$cond==cond1,]$fc),
  #         data.frame(
  #           gene=pvalcalc[pvalcalc$cond==cond2,]$gene,
  #           p2=-log10(pvalcalc[pvalcalc$cond==cond2,]$p),
  #           fc2=pvalcalc[pvalcalc$cond==cond2,]$fc)
  #       )
  #       compared_cond <- merge(compared_cond, genedesc,all.x=TRUE)
  #       compared_cond$genedesc[is.na(compared_cond$genedesc)] <- "Other"
  # 
  #       fc_range <- range(c(compared_cond$fc1, compared_cond$fc2))                
  #       p_range <- range(c(compared_cond$p, compared_cond$p))                
  #       
  #       ggplot(compared_cond, aes(fc1,fc2, label=gene, color=genedesc)) + geom_point(color="gray") + geom_text() +
  #         xlab(paste("FC",cond1)) + ylab(paste("FC",cond2)) + 
  #         xlim(fc_range[1], fc_range[2]) + ylim(fc_range[1], fc_range[2])
  #       ggsave(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "comparisonFC", cond1,"vs",cond2,".pdf")))
  #       ggplot(compared_cond, aes(p1,p2, label=gene, color=genedesc)) + geom_point(color="gray") + geom_text() +
  #         xlab(paste("-log10 pval",cond1)) + ylab(paste("-log10 pval",cond2)) #+
  #         #xlim(p_range[1], p_range[2]) + ylim(p_range[1], p_range[2])
  #       ggsave(file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "comparisonP", cond1,"vs",cond2,".pdf")))
  #     }
  #   }
  # }
  # 
  
  
  #### Plotting of RGR, per mouse 
  list_of_genes <- list()
  for(the_gene in unique(all_gr_meta$gene)){
    #gene <- "PBANKA_100360"
    # onep <- ggplot(all_gr_meta[all_gr_meta$gene==the_gene,],aes(x=day_int,y=rgr, group=mouse_ref, color=condition_ref)) + 
    #   geom_line()+
    #   ggtitle(the_gene)
    onep <- ggplot(all_gr_meta[all_gr_meta$gene==the_gene,],aes(x=day_int,y=rgr, group=mouse_ref, color=mouse_ref)) + 
      geom_line()+
      ggtitle(the_gene)
    onep
    list_of_genes[[the_gene]] <- onep
  }
  ggpubr::ggexport(plotlist=list_of_genes, 
                   filename = file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "allgenes_rgr_permouse.pdf")), 
                   ncol=1, nrow=1)
  
  
  #### Plotting of RGR, per cond 
  list_of_genes <- list()
  for(the_gene in unique(all_gr_meta_percond$gene)){
    onep <- ggplot(all_gr_meta_percond[all_gr_meta_percond$gene==the_gene,],aes(x=day_int,y=rgr, group=condition_ref, color=condition_ref)) + 
      geom_line()+
      ggtitle(the_gene)
    onep
    list_of_genes[[the_gene]] <- onep
  }
  ggpubr::ggexport(plotlist=list_of_genes, 
                   filename = file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "allgenes_rgr_percond.pdf")), 
                   ncol=1, nrow=1)
  
  
  ###################################### Plotting of raw growth curves ... ish ############################
  
  #### Normalize to list of control genes
  normalized_cnt <- cnt
  if(FALSE){
    the_sf <- colMeans(normalized_cnt[
      rownames(normalized_cnt) %in% control_genes$geneid & 
        !(rownames(normalized_cnt) %in% invalid_control_genes),])
    for(i in 1:ncol(normalized_cnt)){
      normalized_cnt[,i] <- normalized_cnt[,i]/the_sf[i]
    }
  }
  
  
  ## Function to plot one gene
  long_cnt <- melt(as.matrix(normalized_cnt)) ######### 666 gives warning
  colnames(long_cnt) <- c("gene","sampleid","exp")
  forgg <- merge(long_cnt,samplemeta[,c("barseqid","sampleid","day_int","mouse_ref","condition_ref")])
  list_of_genes <- list()
  for(the_gene in unique(forgg$gene)){
    onep <- ggplot(forgg[forgg$gene==the_gene,],aes(x=day_int,y=exp, group=mouse_ref, color=condition_ref)) + 
      geom_line()+
      ggtitle(the_gene)
    onep
    list_of_genes[[the_gene]] <- onep
  }
  ggpubr::ggexport(plotlist=list_of_genes, 
                   filename = file.path("/corgi/otherdataset/ellenbushell/plots_barseq_qc", paste(current_pool, "allgenes_normalized2total.pdf")), 
                   ncol=1, nrow=1)
  
  
  all_timecourses[[current_pool]] <- list(
    all_gr_meta_percond=all_gr_meta_percond,
    all_gr_meta=all_gr_meta,
    forgg=forgg,
    genes=unique(forgg$gene)
    #geneplots_rawish=geneplots_rawish
  )
  
}
#Finish all_growthcurves
saveRDS(all_timecourses, "/corgi/websites/shinybarseq/timecourses.rds")






################################################################################
###################### GR stats ################################################
################################################################################




list_grstats <- list()
for(current_pool in list_pools){
  print("GR stats")
  print(current_pool)
  
  
  control_gene_file <- file.path("/corgi/otherdataset/ellenbushell/barseq_pools",current_pool,"/controls.csv")
  if(file.exists(control_gene_file)){
    print("using control gene list")
    control_genes <- read.csv(control_gene_file)
    control_genes$reason <- "control"
    normalize_by_geneids <- control_genes$geneid
  } else {
    print("using default control genes")
    control_genes <- read.csv("/corgi/otherdataset/ellenbushell/dispensible_controls_for_barseq.csv",sep="\t")
    normalize_by_geneids <- control_genes$geneid#[control_genes$reason=="Dispensable"]
  }
  
  
  samplemeta <- read.csv(file.path("/corgi/otherdataset/ellenbushell/filtered_barseq_counts",paste(sep="",current_pool, ".filtered.samplemeta.csv")), row.names = "X")
  cnt <- read.csv(file.path("/corgi/otherdataset/ellenbushell/filtered_barseq_counts",paste(sep="",current_pool, ".filtered.counts.csv")), row.names = "X")
  
  cnt <- cnt[rowMeans(cnt)>100,]
  print(dim(cnt))
  
  #Ensure days in order
  samplemeta$day_int <- as.integer(samplemeta$day)  
  cnt <- cnt[,order(samplemeta$day_int)]
  samplemeta <- samplemeta[order(samplemeta$day_int),]
  
  #Only keep sensible samples; no inputs
  keep <- !samplemeta$is_input & !is.na(samplemeta$day) & samplemeta$genotype!="ug"
  samplemeta <- samplemeta[keep,]
  cnt <- cnt[,keep]
  print(dim(cnt))
  
  #Normalize each sample based on total count
  for(i in 1:ncol(cnt)){
    cnt[,i] <- (1+cnt[,i])/sum(cnt[,i])
  }
  
  #Normalize each sample based on total count of control genes - output becomes garbage
  # for(i in 1:ncol(cnt)){
  #   cnt[,i] <- (1+cnt[,i])/sum(cnt[rownames(cnt) %in% normalize_by_geneids,i])
  # }
  
  #Give each mouse an ID
  samplemeta$mouse_ref <- sprintf("%s_%s_%s_%s", samplemeta$pool, samplemeta$mouse, samplemeta$genotype, samplemeta$primed)
  samplemeta$mouse_ref[samplemeta$is_input] <- ""
  
  #Give each condition an ID    
  samplemeta$condition_ref <- sprintf("%s_%s_%s", samplemeta$pool, samplemeta$genotype, samplemeta$primed)
  
  
  #### Compute GR for each strain
  all_gr <- NULL
  alldays <- unique(samplemeta$day_int)
  for(cur_mr in unique(samplemeta$mouse_ref)){
    
    #For each day
    for(curday in alldays[-1]){
      
      #Figure out which samples to compare. Ignore cases of missing samples        
      this_samp <- which(samplemeta$mouse_ref==cur_mr & samplemeta$day_int==curday)
      prev_samp <- which(samplemeta$mouse_ref==cur_mr & samplemeta$day_int==curday-1)
      if(length(this_samp)==1 & length(prev_samp)==1){
        
        all_gr <- rbind(
          all_gr,
          data.frame(
            mouse_ref=cur_mr,
            condition_ref=samplemeta$condition_ref[this_samp],
            gene=rownames(cnt),
            day_int=curday,
            cnt_this=cnt[,this_samp],
            cnt_prev=cnt[,prev_samp]
          )
        )
      }
    }  
  }
  
  ## Log growth rate is on a more sensible scale
  all_gr$log_gr <- log2(all_gr$cnt_this / all_gr$cnt_prev)
  ## Log minimum count is later correlated to dispersion
  all_gr$depth <- log(pmin(all_gr$cnt_this, all_gr$cnt_prev))
  
  ###### Fit dispersion for all GRs
  
  #Sort from high to low. Needed to later ensure variance increases down depth
  all_gr <- all_gr[order(all_gr$depth, decreasing = TRUE),]
  
  #Fit local GR and compensate    
  fitted_lm <- lm(log_gr ~ depth, all_gr)
  all_gr$adjusted_log_gr <- all_gr$log_gr - as.double(fitted_lm$fitted.values)  ## NA values in fitted model
  
  #Plot the fitted mean, which is likely the bias from the +1 pseudocount
  if(TRUE){
    png(file.path(dir_fancyplot, paste(current_pool, "fitted_mean.png")))
    #pdf(file.path(dir_fancyplot, paste(current_pool, "fitted_mean.pdf")))
    plot(all_gr$log_gr, all_gr$depth, pch=19)
    lines(c(0,0),c(-12,0), col="red")
    points(as.double(fitted_lm$fitted.values), all_gr$depth, col="green", pch=19)
    dev.off()
  }
  
  #Each day has different dispersion, so fit individually
  all_gr$local_sd <- NA
  for(one_day in unique(all_gr$day_int)){ ######### 666 crash
    for_gr <- all_gr$day_int==one_day
    #for_gr <- all_gr$day_int>0
    one_gr <- all_gr[for_gr,]
    
    #Fit local GR and compensate    
    #fitted_lm <- lm(log_gr ~ depth, one_gr)
    #all_gr$adjusted_log_gr[for_gr] <- all_gr$log_gr[for_gr] - as.double(fitted_lm$fitted.values)
    
    #Fit local sd
    one_gr$local_sd <- NA
    for(i in 1:nrow(one_gr)){
      use_bw <- 2
      one_gr$local_sd[i] <- sd(one_gr$adjusted_log_gr[
        one_gr$depth > one_gr$depth[i]-use_bw & 
          one_gr$depth < one_gr$depth[i]+use_bw])  #can end up NA
    }
    
    #Assumption: SD continuously increases for lower depths
    #Work down in order and propagate the SD
    cur_sd <- one_gr$depth[1] 
    for(i in 1:nrow(one_gr)){
      if(!is.na(one_gr$local_sd[i]) & one_gr$local_sd[i] > cur_sd){  ##### if NA value, problem!  ############################ crash here!
        cur_sd <- one_gr$local_sd[i]
      }
      one_gr$local_sd[i] <- cur_sd
    }
    
    #Transfer back to "big" GR table
    all_gr$local_sd[for_gr] <- one_gr$local_sd
  }
  
  #Plot adjust gr and fitted sd    
  if(TRUE){
    png(file.path(dir_fancyplot, paste(current_pool, "fitted_dispersion.png")))
    #pdf(file.path(dir_fancyplot, paste(current_pool, "fitted_dispersion.pdf")))
    plot(all_gr$adjusted_log_gr, all_gr$depth, pch=19)
    lines(c(0,0),c(-12,0))
    points(all_gr$local_sd, all_gr$depth, pch=19, col="green")
    dev.off()
    #ggplot(all_gr, aes(local_sd, depth, group=day_int, color=day_int)) + geom_line()
  }
  
  #Shrink the GR based on dispersion
  all_gr$adjusted_log_gr <- all_gr$adjusted_log_gr / all_gr$local_sd
  ggplot(all_gr, aes(adjusted_log_gr, depth, group=day_int, color=day_int)) + geom_point()
  ggsave(file.path(dir_fancyplot, paste(current_pool, "shrunk_GR.png")), dpi=100)
  
  #Compute GR per mouse; weighted by expected dispersion
  all_gr_permouse <- NULL
  for(mr in unique(all_gr$mouse_ref)){
    for(one_gene in unique(all_gr$gene)){
      one_gr <- all_gr[all_gr$mouse_ref==mr & all_gr$gene==one_gene,]
      #https://en.wikipedia.org/wiki/Inverse-variance_weighting
      use_weight <- 1/(one_gr$local_sd**2)
      weighted_avg <- sum(one_gr$adjusted_log_gr*use_weight) / sum(use_weight)
      weighted_var <- 1/sum(use_weight) #assuming each var[gr]=1
      
      all_gr_permouse <- rbind(all_gr_permouse, data.frame(
        mouse_ref=mr,
        condition_ref=one_gr$condition_ref[1],
        gene=one_gene,
        gr_mean=weighted_avg,
        gr_var=weighted_var
      ))
    }
  }
  
  #Compute GR per condition
  all_gr_percond <- NULL
  for(cur_cond in unique(all_gr_permouse$condition_ref)){
    for(one_gene in unique(all_gr_permouse$gene)){
      one_gr <- all_gr_permouse[all_gr_permouse$condition_ref==cur_cond & all_gr_permouse$gene==one_gene,]
      all_gr_percond <- rbind(all_gr_percond, data.frame(
        cond=cur_cond, #one_gr$condition_ref[1],
        gene=one_gene,
        gr_mean=mean(one_gr$gr_mean),
        gr_var=sum(one_gr$gr_var)/(nrow(one_gr)**2)
      ))
    }
  }
  
  #mouse_ref & condition_ref, one value each
  # nrow((all_gr_permouse[all_gr_permouse$gene=="PBANKA_051490",c("mouse_ref","condition_ref","gene")])) #and they are very different
  # nrow(unique(all_gr_percond[all_gr_percond$gene=="PBANKA_051490",c("cond","gene")])) #and they are very different
  # nrow((all_gr_percond[all_gr_percond$gene=="PBANKA_051490",c("cond","gene")])) #and they are very different
  # foo <- (all_gr_percond[all_gr_percond$gene=="PBANKA_051490",c("cond","gene")])
  # foo
  
  #Turn it into RGR
  for(cur_cond in unique(all_gr_percond$cond)){
    
    #Find mean and var of control genes    
    control_gr <- all_gr_percond[all_gr_percond$cond %in% cur_cond & all_gr_percond$gene %in% normalize_by_geneids,]
    control_mean <- mean(control_gr$gr_mean)
    control_var <- sum(control_gr$gr_var)/(nrow(control_gr)**2)
    print(paste(cur_cond, "to RGR",control_mean,control_var))
    
    #Add to all genes
    all_gr_percond$gr_mean[all_gr_percond$cond %in% cur_cond] <- all_gr_percond$gr_mean[all_gr_percond$cond %in% cur_cond] - control_mean
    all_gr_percond$gr_var[all_gr_percond$cond %in% cur_cond] <- all_gr_percond$gr_var[all_gr_percond$cond %in% cur_cond] + control_var
  }

  #If multiple pools, need to merge into one mean & var
  all_gr_percond <- merge(all_gr_percond,unique(data.frame(
    cond=samplemeta$condition_ref,
    gp=paste(samplemeta$genotype,samplemeta$primed,sep="_")
  )))
  all_gr_percond$cond <- all_gr_percond$gp
  all_gr_percond <- sqldf("select cond, gene, avg(gr_mean) as gr_mean, avg(gr_var)/sqrt(count(gr_var)) as gr_var from all_gr_percond group by cond,gene")
  head(all_gr_percond)
  
  
  #Estimate p-values, assuming a normal distribution
  all_gr_percond$p <- pnorm(all_gr_percond$gr_mean,mean=0,sd=sqrt(all_gr_percond$gr_var))
  all_gr_percond$p[all_gr_percond$p > 0.5] <- 1 - all_gr_percond$p[all_gr_percond$p > 0.5]
  all_gr_percond$log_p <- -log10(all_gr_percond$p)
  
  #duplicates in here all_gr_percond ################# this would cause duplicates later!
  all_gr_percond[all_gr_percond$gene=="PBANKA_051490",] #and they are very different
  
  #### Compare all conditions -- scatter and volcano plots
  list_cond_compare <- data.frame(
    cond1=c("BL6_P",    "BL6_NP"),
    cond2=c("RAG1KO_P", "RAG1KO_NP")
  )
  list_scatter <- list()
  for(curcondi in 1:nrow(list_cond_compare)){
    cond1 <- list_cond_compare$cond1[curcondi]
    cond2 <- list_cond_compare$cond2[curcondi]
    samp1 <- str_ends(all_gr_percond$cond, cond1)  #### here we get duplicates!!! or? ends is dangerous. should loop over pools
    samp2 <- str_ends(all_gr_percond$cond, cond2)   
    
    if(sum(samp1)>0 & sum(samp2)>0){
      
      toplot <- merge(
        data.frame(
          gene=all_gr_percond[samp1,]$gene,
          var1=all_gr_percond[samp1,]$gr_var,
          p1=all_gr_percond[samp1,]$log_p,
          fc1=all_gr_percond[samp1,]$gr_mean),
        
        data.frame(
          gene=all_gr_percond[samp2,]$gene,
          var2=all_gr_percond[samp2,]$gr_var,
          p2=all_gr_percond[samp2,]$log_p,
          fc2=all_gr_percond[samp2,]$gr_mean)
      )
      toplot <- merge(toplot, genedesc,all.x=TRUE)
      toplot$genedesc[is.na(toplot$genedesc)] <- "Other"
      
      toplot$normalizedby <- "no"
      toplot$normalizedby[toplot$gene %in% normalize_by_geneids] <- "yes"
      
      fc_range <- range(c(toplot$fc1, toplot$fc2))
      p_range <- range(c(toplot$p1, toplot$p))
      
      ggplot(toplot, aes(fc1,fc2, label=gene, color=genedesc)) + geom_point(color="gray") + geom_text(size=1) +
        xlab(paste("FC",cond1)) + ylab(paste("FC",cond2)) +
        xlim(fc_range[1], fc_range[2]) + ylim(fc_range[1], fc_range[2])
      ggsave(file.path(dir_fancyplot, paste(current_pool, "comparisonFC", cond1,"vs",cond2,".pdf")))
      
      fig <- plot_ly(type = 'scatter', mode = 'markers') 
      fig <- fig %>%
        add_trace(
          x = toplot$fc1,    #not found!
          y = toplot$fc2, 
          text= toplot$gene,
          hoverinfo = 'text',
          color = toplot$normalizedby,
          showlegend = F
        ) %>%
        layout(
          xaxis = list(title = paste("FC",cond1)), 
          yaxis = list(title = paste("FC",cond2)))
      #fig
      #      htmlwidgets::saveWidget(fig, file.path(dir_fancyplot, paste(current_pool, "comparisonFC", cond1,"vs",cond2,".html")), selfcontained = F, libdir = "lib")
      
      
      if(FALSE){
        ggplot(toplot, aes(p1,p2, label=gene, color=genedesc)) + geom_point(color="gray") + geom_text(size=1) +
          xlab(paste("-log10 pval",cond1)) + ylab(paste("-log10 pval",cond2)) #+
        #xlim(p_range[1], p_range[2]) + ylim(p_range[1], p_range[2])
        ggsave(file.path(dir_fancyplot, paste(current_pool, "comparisonP", cond1,"vs",cond2,".pdf")))
      }

      
      #Estimate p-value of difference, assuming a normal distribution
      toplot$diff_fc <- toplot$fc1 - toplot$fc2
      toplot$diff_var <- toplot$var1 + toplot$var2
      toplot$diff_p <- pnorm(toplot$diff_fc, mean=0, sd=sqrt(toplot$diff_var))
      toplot$diff_p[toplot$p > 0.5] <- 1 - toplot$p[toplot$p > 0.5]
      toplot$diff_log_p <- -log10(toplot$diff_p)
     
      list_scatter[[paste(cond1,"/",cond2)]] <- toplot
    }
  }
  
  
  ### volcano plots for each condition
  print("make volcanoes")
  list_volcano <- list()
  for(cond1 in unique(all_gr_percond$cond)){
    toplot <- data.frame(
      gene=all_gr_percond[all_gr_percond$cond==cond1,]$gene,
      p1=all_gr_percond[all_gr_percond$cond==cond1,]$log_p,
      fc1=all_gr_percond[all_gr_percond$cond==cond1,]$gr_mean
    )
    toplot <- merge(toplot, genedesc,all.x=TRUE)
    toplot$genedesc[is.na(toplot$genedesc)] <- "Other"
    
    toplot$normalizedby <- "no"
    toplot$normalizedby[toplot$gene %in% normalize_by_geneids] <- "yes"
    
    ggplot(toplot, aes(fc1, p1, label=gene, color=genedesc)) + geom_point(color="gray") + geom_text(size=1) +
      xlab(paste("FC",cond1)) + ylab(paste("-log10 pval",cond1))
    ggsave(file.path(dir_fancyplot, paste(current_pool, "volcano", cond1,".pdf")))
    
    
    fig <- plot_ly(type = 'scatter', mode = 'markers') 
    fig <- fig %>%
      add_trace(
        x = toplot$fc1, 
        y = toplot$p1, 
        text= toplot$gene,
        hoverinfo = 'text',
        color = toplot$normalizedby,
        #marker = list(color='green'),
        showlegend = F
      ) %>%
      layout(
        xaxis = list(title = paste("FC",cond1)), 
        yaxis = list(title = paste("-log10 pval",cond1)))
    fig
    #    htmlwidgets::saveWidget(fig, file.path(dir_fancyplot, paste(current_pool, "volcano", cond1,".html")), selfcontained = F, libdir = "lib")
    
    list_volcano[[cond1]] <- toplot
    
  }
  
  
  
  # ### Scatter plot comparing two conditions: BL6 / RAG1KO
  # print("make scatter plot")
  # list_scatter2 <- list()
  # if(TRUE){
  #   samp1 <- str_ends(all_gr_percond$cond, "BL6_P")
  #   samp2 <- str_ends(all_gr_percond$cond, "RAG1KO_P")
  #   
  #   if(sum(samp1)>0 & sum(samp2)>0){
  #     
  #     toplot <- merge(
  #       data.frame(
  #         gene=all_gr_percond[samp1,]$gene,
  #         var1=all_gr_percond[samp1,]$gr_var,
  #         fc1=all_gr_percond[samp1,]$gr_mean),
  #       
  #       data.frame(
  #         gene=all_gr_percond[samp2,]$gene,
  #         var2=all_gr_percond[samp2,]$gr_var,
  #         fc2=all_gr_percond[samp2,]$gr_mean)
  #     )
  #     toplot <- merge(toplot, genedesc,all.x=TRUE)
  #     toplot$genedesc[is.na(toplot$genedesc)] <- "Other"
  #     
  #     toplot$delta_fc <- toplot$fc1 - toplot$fc2
  #     toplot$total_var <- toplot$var1 + toplot$var2
  #     
  #     #Estimate p-values, assuming a normal distribution
  #     toplot$p <- pnorm(toplot$delta_fc,mean=0,sd=sqrt(toplot$total_var))
  #     toplot$p[toplot$p > 0.5] <- 1 - toplot$p[toplot$p > 0.5]
  #     toplot$log_p <- -log10(toplot$p)
  #     
  #     toplot$normalizedby <- "no"
  #     toplot$normalizedby[toplot$gene %in% normalize_by_geneids] <- "yes"
  #     
  #     ggplot(toplot, aes(delta_fc,log_p, label=gene, color=genedesc)) + geom_point(color="gray") + geom_text(size=1) +
  #       xlab(paste("FC BL6 P - FC RAG1KO P")) + ylab(paste("log10 -pval")) 
  #     ggsave(file.path(dir_fancyplot, paste(current_pool, "comparisonDELTA P.pdf")))
  #     
  #     
  #     fig <- plot_ly(type = 'scatter', mode = 'markers') 
  #     fig <- fig %>%
  #       add_trace(
  #         x = toplot$delta_fc, 
  #         y = toplot$log_p, 
  #         text= paste(toplot$gene, toplot$genedesc),
  #         hoverinfo = 'text',
  #         color=toplot$normalizedby,
  #         #marker = list(color='green'),
  #         showlegend = F
  #       ) %>%
  #       layout(
  #         xaxis = list(title = paste("FC BL6 P - FC RAG1KO P")), 
  #         yaxis = list(title = paste("log10 -pval")))
  #     fig
  #     #      htmlwidgets::saveWidget(fig, file.path(dir_fancyplot, paste(current_pool, "comparisonDELTA P.html")), selfcontained = F, libdir = "lib")
  #     
  #     
  #     list_scatter2[["BL6_P vs RAG1KO_p"]] <- toplot
  #     
  #   }
  #}
  
  list_grstats[[current_pool]] <- list(
    volcano=list_volcano,
    scatterplot=list_scatter
    #scatterplot2=list_scatter2
  )

}
saveRDS(list_grstats, "/corgi/websites/shinybarseq/grstats.rds")
#Finish GRstats


if(FALSE){
  list_grstats <- readRDS("/corgi/websites/shinybarseq/grstats.rds")
  #list_grstats$EB_barseq_slowpool_1$
}
