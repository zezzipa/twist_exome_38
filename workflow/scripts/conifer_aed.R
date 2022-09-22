#CREATE AED
if (file.size(snakemake@input[["txtfile"]]) == 0) {
  df <- "Conifer did not find any structural variants."
  write.table(df, file = snakemake@output[["aedfile"]], row.names = FALSE, quote = FALSE, sep = "\t")
} else {
  df <- read.csv(snakemake@input[["txtfile"]], header=FALSE, sep = "\t", stringsAsFactors=FALSE)

  df$V1 <- sub('.+:(.+)', '\\1', df$V1)

  write.table(df, file = snakemake@input[["txtfile"]], col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

  cNames <- c("name","chromosome","start","end","type")
  names(df) <- cNames
  keep <- c("chromosome","start","end","name","type","type")
  df = df[keep]

  df$type[df$type=="dup"]<-"copynumber/gain"
  df$type[df$type=="del"]<-"copynumber/loss"
  df$type.1[df$type.1=="dup"]<-"rgb(0,0,255)"
  df$type.1[df$type.1=="del"]<-"rgb(255,0,0)"

  df$new <- NA #blank column
  NewNames <- c("bio:sequence(aed:String)","bio:start(aed:Integer)","bio:end(aed:Integer)","aed:name(aed:String)","aed:category(aed:String)","style:color(aed:Color)","aed:value(aed:String)")
  names(df) <- NewNames
  df <- df[, c("bio:sequence(aed:String)","bio:start(aed:Integer)","bio:end(aed:Integer)","aed:name(aed:String)","aed:value(aed:String)","aed:category(aed:String)","style:color(aed:Color)")]

  #HEADER
  header2 <- data.frame("","","","affx:ucscGenomeVersion(aed:String)","hg19","","", stringsAsFactors=FALSE)
  names(header2) <- c("bio:sequence(aed:String)","bio:start(aed:Integer)","bio:end(aed:Integer)","aed:name(aed:String)","aed:value(aed:String)","aed:category(aed:String)","style:color(aed:Color)")
  df <- rbind(header2, df)
  header1 <- data.frame("","","","namespace:affx(aed:URI)","http://affymetrix.com/ontology/","","", stringsAsFactors=FALSE)
  names(header1) <- c("bio:sequence(aed:String)","bio:start(aed:Integer)","bio:end(aed:Integer)","aed:name(aed:String)","aed:value(aed:String)","aed:category(aed:String)","style:color(aed:Color)")
  df <- rbind(header1, df)

  write.table(df, file = snakemake@output[["aedfile"]], row.names = FALSE, quote = FALSE, sep = "\t")
}
