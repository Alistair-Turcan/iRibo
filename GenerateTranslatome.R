time1 = Sys.time()

# Function to parse command line arguments
parseArguments <- function(args) {
  arguments <- list()
  
  for(i in 2:length(args)) {
    arg <- args[[i]]
    
    if(substr(arg, 1, 2) == "--") {
      delimiterPos <- regexpr("=", arg)
      if(delimiterPos != -1) {
        key <- substr(arg, 3, delimiterPos - 1)
        value <- substr(arg, delimiterPos + 1, nchar(arg))
        arguments[[key]] <- value
		print(key)
      } else {
        key <- substr(arg, 3, nchar(arg))
        arguments[[key]] <- ""
      }
    }
  }
  
  return(arguments)
}

# Function to retrieve argument value with a default or stop execution if mandatory and not specified
getArg <- function(arguments, arg, defaultValue, mandatory = FALSE) {
  if(arg %in% names(arguments)) {
    return(arguments[[arg]])
  } else {
    if(mandatory) {
      stop(paste("Argument", arg, "is required but not specified."))
    } else {
      return(defaultValue)
    }
  }
}



args <- commandArgs(trailingOnly = FALSE)
parsedArgs <- parseArguments(args)

calls_file <- getArg(parsedArgs, "TranslationCalls", "", TRUE)
orfs_file <- getArg(parsedArgs, "CandidateORFs", "", TRUE)

threads <- as.integer(getArg(parsedArgs, "Threads", "1", FALSE))

exclude_chr_input <- getArg(parsedArgs, "ExcludeCHR", "", FALSE)
exclude_chr <- strsplit(exclude_chr_input, ",")[[1]]

exclude_overlap_gene <- tolower(getArg(parsedArgs, "ExcludeOverlapGene", "", FALSE))

output_dir <- getArg(parsedArgs, "Output", "", FALSE)

if(output_dir!=""){
	dir.create(output_dir, showWarnings = FALSE)

	output_dir = paste(output_dir, "/", sep = "")
}

#Load in calls file from TranslationCalls, and original list of orfs from GetCandidateORFs
calls<-read.csv(calls_file,sep=" ",stringsAsFactors=F)
orfs<-read.csv(orfs_file, sep=" ",stringsAsFactors=F)

FDR  <- as.numeric(getArg(parsedArgs, "FDR", "0.05", FALSE))

library(scales)
library(parallel)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
	install.packages("codetools")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

if (!requireNamespace("future.apply", quietly = TRUE)) {
  install.packages("future.apply")
}
library("future.apply")
plan(multicore, workers = threads) #Enables multithreading

print(Sys.time() - time1)
time1 = Sys.time()


# Define a function to calculate the binom.test for a single index i
#For scrambled ORFs
calc_binom_false <- function(i) {
  temp_call <- calls[i,]
  scrambled_sum = temp_call$scrambled_sum
  #scrambled_sum =  temp_call$scrambled0 + temp_call$scrambled1 + temp_call$scrambled2
  if ( scrambled_sum> 0) {
		return(binom.test(temp_call$scrambled0, scrambled_sum, 1/3, alt = "g")[[3]])
 } else {
    return(1)
  }
}

# Export the calls object to the worker R processes
clusterExport(cl = makeCluster(threads), varlist = c("calls"))

# Use future_lapply to apply the calc_binom_false function to each index in parallel
results <- future_lapply(1:length(orfs[,1]), function(i) calc_binom_false(i))

print(Sys.time() - time1)
time1 = Sys.time()

# Combine the results into a single vector
scram_bin <- unlist(results)

print(Sys.time() - time1)
time1 = Sys.time()

# Define a function to calculate the binom.test for a single index i
#For original, unscrambled ORFs
calc_binom_true <- function(i) {
  temp_call <- calls[i,]
  frame_sum=temp_call$frame_sum
  #frame_sum=temp_call$frame0 + temp_call$frame1 + temp_call$frame2
  if ( temp_call$frame_sum > 0) {

		return(binom.test(temp_call$frame0, frame_sum, 1/3, alt = "g")[[3]])
	

 } else {
    return(1)
  }
}


# Use future_lapply to apply the calc_binom_false function to each index in parallel
results <- future_lapply(1:length(orfs[,1]), function(i) calc_binom_true(i))

# Combine the results into a single vector
ribo_bin <- unlist(results)
print(Sys.time() - time1)
time1 = Sys.time()

#Identify indicies of canonical and noncanonical ORFs
#overlaps<-read.csv("orf_overlap",sep=" ") #a file giving info on what annotation each ORF overlaps
#ORFs are candidates if they do not overlap an annotated gene on the same strand
#all_index<-which(orfs$splice_gene=="X" & orfs$contig<16 & orfs$is_gene!="YDL185W") #ORFs overlapping spliced genes are excluded. "YDL185W" contains entein and is excluded. We also exclude ORFs on contig 16, which is mitochondria
#all_index<-which(orfs$chr<16 & orfs$gene_id!="YDL185W") #ORFs overlapping spliced genes are excluded. "YDL185W" contains entein and is excluded. We also exclude ORFs on contig 16, which is mitochondria

#candidate_index<-all_index[which((overlaps$sense_ver+overlaps$sense_unchar+overlaps$sense_te+overlaps$sense_blocked)[all_index]==0)]

canonical_index <- which(orfs$gene_id != "")
#canonical_index <- which(orfs$orf_class %in% c("Verified", "Uncharacterized", "transposable_element_gene"))

noncanonical_index <- which(orfs$gene_id == "")
if(exclude_overlap_gene == "true"){
	print("Excluding Overlapping nORFs")
	noncanonical_index <- intersect(noncanonical_index, which(orfs$CDS_intersect == "X"))
}
print(length(noncanonical_index))
noncanonical_index <- intersect(noncanonical_index, which(!(orfs$chr_str %in% exclude_chr)))
print(length(noncanonical_index))

#noncanonical_index <- candidate_index[which(!(orfs$orf_class[candidate_index] %in% c("Verified", "Uncharacterized", "transposable_element_gene")))]

print(length(noncanonical_index))
print(length(canonical_index))



print(Sys.time() - time1)
time1 = Sys.time()

num_hits_noncanonical<-array()
scrambled_hits_noncanonical<-array()
num_hits_canonical<-array()
scrambled_hits_canonical<-array()

#number of true hits and scrambled (negative control) hits at range of pval thresholds, used to calculate FDR

# Define a function to calculate the hits and scrambled hits at a given pval threshold
calc_hits <- function(i, ribo_bin, scram_bin, index) {
  num_hits <- length(which(ribo_bin[index] < i/10000))
  scrambled_hits <- length(which(scram_bin[index] < i/10000))
  return(list(num_hits = num_hits, scrambled_hits = scrambled_hits))
}

# Use future_lapply to apply the calc_hits function to each pval threshold in parallel
results <- future_lapply(1:2000, function(i) {
  num_noncanonical <- calc_hits(i, ribo_bin, scram_bin, noncanonical_index)
  num_canonical <- calc_hits(i, ribo_bin, scram_bin, canonical_index)
  return(list(num_noncanonical = num_noncanonical, num_canonical = num_canonical))
})

# Combine the results into arrays
for(i in 1:2000) {
  temp_results = results[[i]]
  num_hits_noncanonical[i] <- temp_results$num_noncanonical$num_hits
  scrambled_hits_noncanonical[i] <- temp_results$num_noncanonical$scrambled_hits
  num_hits_canonical[i] <- temp_results$num_canonical$num_hits
  scrambled_hits_canonical[i] <- temp_results$num_canonical$scrambled_hits
}

print(Sys.time() - time1)
time1 = Sys.time()

df_hits<-data.frame(pvals=(1:2000)/10000,hits=c(num_hits_noncanonical,scrambled_hits_noncanonical),hit_type=c(rep("Actual",length(num_hits_noncanonical)),rep("Scrambled",length(scrambled_hits_noncanonical))))
df_hits_canonical<-data.frame(pvals=(1:2000)/10000,hits=c(num_hits_canonical,scrambled_hits_canonical),hit_type=c(rep("Actual",length(num_hits_canonical)),rep("Scrambled",length(scrambled_hits_canonical))))

print(Sys.time() - time1)
time1 = Sys.time()

##### FIGURE SIZE
one.c <- 85 #single column
one.5c <- 144 #1.5 column
two.c <- 174 #full width

##### TEXT SIZE
titles <- 8*2
txt <- 8*2
lbls <- 9*2

my_p = which.min(abs(scrambled_hits_noncanonical/num_hits_noncanonical-FDR))/10000 - .0001
my_p_canon = which.min(abs(scrambled_hits_canonical/num_hits_canonical-FDR))/10000 - .0001

orfs_found = num_hits_noncanonical[my_p*10000]
orfs_found_canon = num_hits_canonical[my_p_canon*10000]

my_title = paste("P-value: ", my_p)
my_title = paste(my_title, " Orfs Found: ")
my_title = paste(my_title, orfs_found)

f1c<- ggplot(df_hits,aes(x=pvals,y=hits,linetype=hit_type))+
  geom_line(size=1,col='blue')+
  scale_color_manual(values=c("#7C4DFF","#9C27B0"))+
  labs(x="P-value threshold",y="Translated nORFs found")+
  geom_vline(xintercept=which.min(abs(scrambled_hits_noncanonical/num_hits_noncanonical-FDR))/10000, linetype="dashed", color="black", size=1.0)+
  scale_y_continuous() +
  theme_linedraw()+
  ggtitle(my_title) +
  theme(plot.title = element_text(size = titles, face = "bold"),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position=c(.8,.7),
        legend.background = element_rect(fill = 'transparent'))
f1c

print(Sys.time() - time1)
time1 = Sys.time()

plot_output = paste(output_dir, "nORF Discovery.png", sep = "")
png(plot_output)
plot(f1c)
dev.off()

noncanonical_trans<-noncanonical_index[which(ribo_bin[noncanonical_index]<.008)]

print(Sys.time() - time1)
time1 = Sys.time()

#Print the indices to a file
noncan_indices <- noncanonical_index[which(ribo_bin[noncanonical_index] < my_p)] - 1
write.table(noncan_indices, "orf_indices", col.names=FALSE, row.names=FALSE, quote=FALSE)

canon_indices <- canonical_index[which(ribo_bin[canonical_index] < my_p_canon)] - 1
write.table(canon_indices, "orf_indices_can", col.names=FALSE, row.names=FALSE, quote=FALSE)

# Write each ORF that is indicated in noncan_indices to a file
translated_nORFs <- orfs[noncan_indices, ]

# Extract 'in-frame reads' values corresponding to the noncan_indices
in_frame_reads <- sapply(noncan_indices + 1, function(i) calls[i,]$reads0)

# Add 'in-frame reads' column to the dataframe
translated_nORFs$in_frame_reads <- in_frame_reads

##DO THE SAME FOR cORFs

# Add 'Expression Level' column to the dataframe
translated_nORFs$Expression_Level <- translated_nORFs$in_frame_reads / translated_nORFs$orf_length

nORF_output = paste(output_dir, "translated_nORFs.csv", sep = "")
write.csv(translated_nORFs, nORF_output, row.names=FALSE)

# Write each ORF that is indicated in canon_indices to a file
translated_cORFs <- orfs[canon_indices, ]

# Extract 'in-frame reads' values corresponding to the canon_indices
in_frame_reads <- sapply(canon_indices + 1, function(i) calls[i,]$reads0)

# Add 'in-frame reads' column to the dataframe
translated_cORFs$in_frame_reads <- in_frame_reads

# Add 'Expression Level' column to the dataframe
translated_cORFs$Expression_Level <- translated_cORFs$in_frame_reads / translated_cORFs$orf_length

cORF_output = paste(output_dir, "translated_cORFs.csv", sep = "")
write.csv(translated_cORFs, cORF_output, row.names=FALSE)
