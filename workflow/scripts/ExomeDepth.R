library(ExomeDepth)
library(rlang)

{
data(genes.hg19)
data(exons.hg19)
data(Conrad.hg19)

# bedfile maleRefCount.mat bamfiles
#References CPRefCount.mat
load(args[1])

#Probes as bed.frame
probes.hg19 <- read.csv(args[0], header=F, sep = "\t")

# Create counts dataframe for all BAMs
my.bam <- args[-(1:2)]
my.counts <- getBamCounts(bam.files = my.bam,
                          include.chr = F)

ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')

#If bam-files write chr1 etc
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$space),
                                   pattern = 'chr',
                                   replacement = '') ##remove "chr" letters

#Write and move count table
seqid=args[2]
write.table(ExomeCount.dafr, paste("SV/", "ExomeCount_", seqid, ".txt", sep=""), sep = '\t')

### prepare the main matrix of read count data
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = '.bam$')])
nsamples <- ncol(ExomeCount.mat)

### start looping over each sample
for (i in 1:nsamples) {
   #### Create the aggregate reference set for this sample
   my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
                                      reference.counts = CPRefCount.mat,
                                      bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                      n.bins.reduced = 10000)
   my.reference.selected <- apply(X = CPRefCount.mat[, my.choice$reference.choice, drop = FALSE],
                                  MAR = 1,
                                  FUN = sum)

   message('Now creating the ExomeDepth object')
   all.exons <- new('ExomeDepth',
                    test = ExomeCount.mat[,i],
                    reference = my.reference.selected,
                    formula = 'cbind(test, reference) ~ 1')


   ################ Now call the CNVs
   all.exons <- CallCNVs(x = all.exons,
                         transition.probability = 10^-4,
                         chromosome = ExomeCount.dafr$chromosome,
                         start = ExomeCount.dafr$start,
                         end = ExomeCount.dafr$end,
                         name = ExomeCount.dafr$names)

##################### Now annotate the ExomeDepth object
#Exons
   exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19$chromosome,
                                                IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
                                                names = exons.hg19$name)
   all.exons <- AnnotateExtra(x = all.exons,
                              reference.annotation = exons.hg19.GRanges,
                              min.overlap = 0.0001,
                              column.name = 'exons')

#Genes
   genes.hg19.GRanges <- GenomicRanges::GRanges(seqnames = genes.hg19$chromosome,
                                                IRanges::IRanges(start=genes.hg19$start,end=genes.hg19$end),
                                                names = genes.hg19$name)
   all.exons <- AnnotateExtra(x = all.exons,
                              reference.annotation = genes.hg19.GRanges,
                              min.overlap = 0.0001,
                              column.name = 'gene')

## Output files
  #NexusSV file

   df = all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),]
   nexus <- c("id", "type")
   nexus_df = df[nexus]
   nexus_df$type[nexus_df$type=="duplication"]<-"CN Gain"
   nexus_df$type[nexus_df$type=="deletion"]<-"CN Loss"

   names(nexus_df) <- c("Chromosome Region", "Event")

   bamfile <- basename(my.bam[i])
   output.file <- paste('SV/ExomeDepth_', bamfile, '_SV.txt', sep = '')
   output.file <- gsub(x = output.file, ".bam", "")
   write.table(nexus_df, file = output.file, row.names = FALSE, quote = FALSE, sep = "\t")

   #Txt file
    bamfile <- basename(my.bam[i])
    output.file <- paste('SV/ExomeDepth_', bamfile, '.txt', sep = '')
    output.file <- gsub(x = output.file, ".bam", "")
    write.table(file = output.file, x = all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),],
                row.names = FALSE, sep = "\t")
  #AED file
    df <- data.frame(all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),], stringsAsFactors = FALSE)
    keep <- c("chromosome", "start", "end", "gene", "nexons", "reads.ratio", "type", "type")
    df = df[keep]
    df$chromosome <- sub("^", "chr", df$chromosome )
    df$type[df$type=="duplication"]<-"copynumber/gain"
    df$type[df$type=="deletion"]<-"copynumber/loss"
    df$type.1[df$type.1=="duplication"]<-"rgb(0,0,255)"
    df$type.1[df$type.1=="deletion"]<-"rgb(255,0,0)"

    df$new <- NA #blank column
    NewNames <- c("bio:sequence(aed:String)", "bio:start(aed:Integer)", "bio:end(aed:Integer)", "aed:name(aed:String)",
                  "bio:markerCount(aed:Integer)", "bio:state(aed:Rational)", "aed:category(aed:String)",
                  "style:color(aed:Color)", "aed:value(aed:String)")
    names(df) <- NewNames
    df <- df[, c("bio:sequence(aed:String)", "bio:start(aed:Integer)", "bio:end(aed:Integer)", "aed:name(aed:String)",
              "aed:value(aed:String)", "bio:markerCount(aed:Integer)", "bio:state(aed:Rational)",
              "aed:category(aed:String)", "style:color(aed:Color)")]

    #HEADER
    header2 <- data.frame("","","","affx:ucscGenomeVersion(aed:String)","hg19","","","","", stringsAsFactors = FALSE)
    names(header2) <- c("bio:sequence(aed:String)","bio:start(aed:Integer)","bio:end(aed:Integer)","aed:name(aed:String)",
                        "aed:value(aed:String)","bio:markerCount(aed:Integer)","bio:state(aed:Rational)",
                        "aed:category(aed:String)","style:color(aed:Color)")
    df <- rbind(header2, df)
    header1 <- data.frame("","","","namespace:affx(aed:URI)","http://affymetrix.com/ontology/","","","","", stringsAsFactors = FALSE)
    names(header1) <- c("bio:sequence(aed:String)","bio:start(aed:Integer)","bio:end(aed:Integer)","aed:name(aed:String)",
                      "aed:value(aed:String)","bio:markerCount(aed:Integer)","bio:state(aed:Rational)",
                      "aed:category(aed:String)","style:color(aed:Color)")
    df <- rbind(header1, df)

    output.file <- gsub(x = output.file, ".txt", ".aed")
    print(output.file)
    write.table(df, file = output.file, row.names = FALSE, quote = FALSE, sep = "\t")
}

}
