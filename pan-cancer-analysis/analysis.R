# Load data.table library
library(data.table)

# Load SNP array data (https://gdc.cancer.gov/about-data/publications/pancanatlas)
data = fread("pan-cancer-analysis/data/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg", sep="\t")

# Load annotation data (https://gdc.cancer.gov/about-data/publications/pancanatlas)
ann = fread("pan-cancer-analysis/data/merged_sample_quality_annotations.tsv", sep="\t")

# Match annotation data
ann_matched = ann[which(ann$aliquot_barcode %in% data$Sample),]

# Which unique cancer types are available?
types = unique(ann_matched$`cancer type`)

# Function to analyse TCGA SNP array data for copy number alterations in specified region
cn_analyse = function(name, chr, pos1, pos2, dir) {
  
  # Initialize overviews
  overview = list()
  summary = NULL
  complete = NULL
  
  # Min and max for given chromosome
  pos_min = min(as.numeric(data$Start[which(data$Chromosome==chr)]))
  pos_max = max(as.numeric(data$Start[which(data$Chromosome==chr)]))
  
  # Iterate through cancer types
  for (type in types) {
    
    # Select samples from given cancer type AND derived from tumour (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes) AND not marked 'do not use'
    samples = ann_matched$aliquot_barcode[which(ann_matched$`cancer type` == type & substr(ann_matched$aliquot_barcode,14,15)%in%paste0("0",1:9) & ann_matched$Do_not_use == FALSE)]
    
    # Initialize results table
    type_overview = NULL
    
    # Iterate through samples
    for (sample_i in 1:length(samples)) {
      
      # Select SNP array data matching with chr and pos1 and pos2
      sample = samples[sample_i]
      sample_data = data[which(data$Sample==sample & data$Chromosome==chr & data$Start<=pos1& data$End>=pos1),]
      
      # In case no segmentation, take mean as CN value
      cn_val = mean(sample_data$Segment_Mean)
      
      # Check if CN value is valid
      if (!is.na(cn_val)) {
        
        # If CN value < -0.3, call as loss
        if (cn_val < (-0.3)) {
          cn_val = "loss"
        }
        # If CN value > 0.3, call as gain
        else if (cn_val > 0.3) {
          cn_val = "gain"
        }
        # If -0.3 <= CN value <= 0.3, call as normal
        else {
          cn_val = "normal"
        }
        
      }
      # Or mark as NA
      else {
        cn_val = NA
      }
      
      # Save results to type_overview
      type_overview = rbind(type_overview, c(sample, cn_val))
      
      # Print progression
      if((sample_i-1)%%25 == 0) {
        print(paste0("Progression: analyzing ", type, " sample ",sample_i,"/",length(samples)))
      }
    }
    
    # Save results
    colnames(type_overview) = c("sample","CNA")
    overview[[type]] = type_overview
    complete = rbind(complete, type_overview)
    
    total = length(which(!is.na(type_overview[,"CNA"])))
    gain = length(which(type_overview[,"CNA"]=="gain"))
    loss = length(which(type_overview[,"CNA"]=="loss"))
    normal = length(which(type_overview[,"CNA"]=="normal"))
    
    summary = rbind(summary,c(type,
                              total,
                              gain/total,
                              loss/total,
                              normal/total
    ))
  }
  
  # Save RDS files
  saveRDS(overview, paste0(dir,"/",name,"-overview.RDS"))
  saveRDS(summary, paste0(dir,"/",name,"-summary.RDS"))
  saveRDS(complete, paste0(dir,"/",name,"-complete.RDS"))
  
}

# Function to read stored cn_analyse data
cn_read = function(name, dir) {
  
  # Read RDS files
  cn_overview <<- readRDS(paste0(dir,"/",name,"-overview.RDS"))
  cn_summary <<- readRDS(paste0(dir,"/",name,"-summary.RDS"))
  cn_complete <<- readRDS(paste0(dir,"/",name,"-complete.RDS"))
  
}

# For all analyses: use GRCh38 locations for genes (https://www.genecards.org/)
# Analysis of TRB gene
cn_analyse(name = "TRB",
           chr = 7,
           pos1 = 142299011,
           pos2 = 142813287,
           dir = "pan-cancer-analysis/res")

# Analysis of TRD gene
cn_analyse(name = "TRD",
           chr = 14,
           pos1 = 22422546,
           pos2 = 22466577,
           dir = "pan-cancer-analysis/res")

# Read and process TCGA cohorts annotation
cohorts = read.csv("pan-cancer-analysis/data/tcga-cohorts.txt", sep="\t", stringsAsFactors = F, check.names = F)
capitalize = function(x) {
  y = strsplit(x, " ")[[1]]
  paste(toupper(substring(y, 1,1)), substring(y, 2),
        sep="", collapse=" ")
}
rownames(cohorts) = cohorts$Cohort

# Function to create frequency plot
create_cn_freq_plot = function(name, main) {
  
  # Set colors
  color_1 = rgb(.25,.25,.25)
  color_2 = rgb(.75,.75,.75)
  graphics.off()
  
  # Read relevant data
  cn_read(name = name,
          dir = "pan-cancer-analysis/res")
  
  # Initiate PNG
  img_file= paste0("pan-cancer-analysis/res/",name,"-CNV.png")
  png(img_file, res=600, width=5000, height=4500)  
  {
    # Set margin
    par(mar=c(5,5,5,5))
    
    # Set limits of axes
    xlim = c(0,10)
    ylim = c(0,33)
    
    # Initiate plot
    plot(type="n",
         x = c(-12,10),
         y = c(0,33),
         axes = F,
         xlab = "",
         ylab = "",
         main = "",
         bty = "l", 
         xaxs = "i", 
         yaxs = "i")
    
    # Draw axes
    xat = seq(from = 0, to = 10, by = 2.5)
    yat = seq(from = 0, to = 33, by = 1)
    segments(xat, 0, xat, 33, col  ="#eeeeee", lwd=1.4, xpd=T)
    segments(xat, 0, xat, -0.5, col  ="#b1b1b1", lwd=1.4, xpd=T)
    text(xat, -0.5, pos = 1, col = "#333333", labels = c("0%","25%","50%","75%","100%"), xpd = T)  
    text(5, -2.5, pos = 1, col = "#333333", labels = "Fraction affected tumours", xpd = T)  
    text(5, 35.5, pos = 1, col = "#333333", labels = main, xpd = T, font=2)  
    
    # Exclude DLBC and LAML from visualisation and order matrix
    freq_matrix = as.matrix(cn_summary[which(cn_summary[,1] != "DLBC" & cn_summary[,1] != "LAML"),])
    freq_matrix = freq_matrix[order(as.numeric(freq_matrix[,5]),decreasing = T),]
    rownames(freq_matrix) = freq_matrix[,1]
    
    # Plot bars per tumour type
    for (i in 1:nrow(freq_matrix)) {
      p1 = as.numeric(freq_matrix[i,3])*10
      p2 = as.numeric(freq_matrix[i,4])*10
      y=i-0.5+2
      rect(0,y-0.4,p1,y+0.4,border="white",col=color_1, lwd=1.4)
      rect(p1,y-0.4,p1+p2,y+0.4,border="white",col=color_2, lwd=1.4)
      l=capitalize(cohorts[rownames(freq_matrix)[i], "Disease Name"])
      l=stringr::str_replace(string = l, pattern = "And ", replacement = "and ")
      text(-0.25,y,adj=1, col="#333333", labels=paste0(l, " (n=", as.numeric(freq_matrix[rownames(freq_matrix)[i], 2]),")"), xpd=T)  
    }
    
    # Plot pan-cancer mean bars
    p1 = mean(as.numeric(freq_matrix[,3]))*10
    p2 = mean(as.numeric(freq_matrix[,4]))*10
    y=1
    rect(0,y-0.4,p1,y+0.4,border="white",col=color_1, lwd=1.4)
    rect(p1,y-0.4,p1+p2,y+0.4,border="white",col=color_2, lwd=1.4)
    text(-0.25,y,adj=1, col="#333333", labels="Pan-cancer (mean)", xpd=T)  
    
    # Print pan-cancer means
    print(p1)
    print(p2)
    print(p1+p2)
    segments(0,0,0,33,xpd=T,col="#B1B1B1",lwd=1.4)
    segments(0,0,10,0,col="#B1B1B1",lwd=1.4,xpd=T)
    y = 1
    rect(5.25,y-.4,5.75,y+.4,col=color_1,border="white",lwd=1.4)
    text(5.75,y,col="#333333", labels = "gain", pos=4)
    
    y = 1
    rect(7.75,y-.4,8.25,y+.4,col=color_2,border="white",lwd=1.4)
    text(8.25,y,col="#333333", labels = "loss", pos=4)
  }
  
  # Finalise and open PNG
  dev.off()
  system(paste0("open ",img_file))
}

# Read TRB data
create_cn_freq_plot("TRB", "") #or use main=substitute(paste(bold("Copy number alterations in "), bolditalic('TRB'))))

# Read TRD data
create_cn_freq_plot("TRD", "") #or use main=substitute(paste(bold("Copy number alterations in "), bolditalic('TRD'))))

# Function to visualise copy number alterations across chromosomes
plot_cn = function(data, chr, pos1, pos2, filename="") {
  
  # Set colors
  color_1 = rgb(.25,.25,.25)
  color_2 = rgb(.75,.75,.75)
  cytoband=read.table("pan-cancer-analysis/data/cytobands.txt",sep="\t",stringsAsFactors = F,check.names = F)
  
  detailed = data
  samples = unique(data$Sample)
  
  chr_pos = c(0,max(cytoband$V3[which(cytoband$V1==paste0("chr",chr))]))
  yscale = length(samples)/10
  plot(chr_pos,
       c(0.5-2*yscale,length(samples)+0.5),
       type="n",
       axes=F,ylab="",xlab="",main="")
  
  y1 = 0.5-yscale
  y2 = 0.5-yscale*1.5
  ym = 0.5-yscale*1.25
  #rect(pos1,0.5-2*yscale,pos2,length(samples)+0.5,col="#333333", border="#333333")
  
  # Iterate through affected tumours to determine order of 
  sums = NULL
  for (y in 1:length(samples)) {
    sample = samples[y]
    sample_data = detailed[which(detailed$Sample==sample),]
    vals = sample_data$Segment_Mean
    vals[which(vals>0.3)] = 1
    vals[which(vals<(-0.3))] = -1
    vals[which(vals>-0.3&vals<(0.3))] = 0
    sums = c(sums,sum(sample_data$Num_Probes*vals*100))
  }
  
  for (y in 1:length(samples)) {
    sample = samples[order(sums)][y]
    sample_data = detailed[which(detailed$Sample==sample),]
    #rect(chr_pos[1],y-0.5,chr_pos[2],y+0.5,col="#DDDDDD",border=NA)
    for (i in 1:nrow(sample_data)) {
      val=sample_data$Segment_Mean[i]
      if (val < (-0.3)) {
        col=color_2
      }
      else if (val > 0.3) {
        col=color_1
      }
      
      else {
        col="white"
      }
      rect(sample_data$Start[i],y-0.5,sample_data$End[i],y+0.5,col=col,border=NA)
    }
  }
  
  # Plot chromosome
  cols=rep("blue",8)
  names(cols) = unique(cytoband$V5)
  cytoband_chr = cytoband[which(cytoband$V1==paste0("chr",chr)),]
  cols["gneg"] = "white"
  cols["gpos25"] = rgb(.75,.75,.75)
  cols["gpos50"] = rgb(.5,.5,.5)
  cols["gpos75"] = rgb(.25,.25,.25)
  cols["gpos100"] = rgb(0,0,0)
  cols["acen"] = "darkred"
  for (i in 1:nrow(cytoband_chr)) { 
    cytoband_data = cytoband_chr[i,]  
    rect(cytoband_data$V2,y1,cytoband_data$V3,y2,col=cols[cytoband_data$V5], border=NA)
  }
  acen = cytoband_chr[which(cytoband_chr$V5=="acen"),]
  rect(min(acen$V2),0.5,max(acen$V3),length(samples)+0.5,col="white", border=NA)
  segments(0,y1,0,y2)
  segments(0,y1,min(acen$V2))
  segments(0,y2,min(acen$V2))
  polygon(c(min(acen$V2),min(acen$V3),min(acen$V3)),c(y1,y1,ym),col="white",border=NA)
  polygon(c(min(acen$V2),min(acen$V3),min(acen$V3)),c(y2,y2,ym),col="white",border=NA)
  segments(min(acen$V2),y1,min(acen$V3),ym)
  segments(min(acen$V2),y2,min(acen$V3),ym)
  segments(chr_pos[2],y1,chr_pos[2],y2)
  segments(chr_pos[2],y1,max(acen$V3))
  segments(chr_pos[2],y2,max(acen$V3))
  polygon(c(max(acen$V2),max(acen$V2),max(acen$V3)),c(ym,y1,y1),col="white",border=NA)
  polygon(c(max(acen$V2),max(acen$V2),max(acen$V3)),c(ym,y2,y2),col="white",border=NA)
  segments(max(acen$V2),ym,max(acen$V3),y1)
  segments(max(acen$V2),ym,max(acen$V3),y2)
  
}

# Visualise TRB copy number alterations
graphics.off()
img_file= "pan-cancer-analysis/res/TRB-plot.png"
png(img_file, res=600, width=5000, height=4500)  
chr = 7
pos1 = 142299011
pos2 = 142813287
detailed_data = data[which(data$Sample%in%cn_complete[which(cn_complete[,"CNA"]!="normal"),"sample"] & data$Chromosome==chr),]
plot_cn(detailed_data,chr,pos1,pos2)
dev.off()