### Install the R phylogenetics package, called "ape"
install.packages("ape", lib="~/genomics-workshop/rlibs", repos="https://cran.r-project.org")
library(ape)

### Add Function for read.vcf HERE!
read.vcf <- function(file, special.char="##", ...) {
    my.search.term=paste0(special.char, ".*")
    all.lines=readLines(file)
    clean.lines=gsub(my.search.term, "",  all.lines)
    clean.lines=gsub("#CHROM", "CHROM", clean.lines)
    read.table(..., text=paste(clean.lines, collapse="\n"))
}

### Then read in the vcf file
my.data=read.vcf("~/genomics-workshop/SNPs/combined.vcf", header=TRUE, stringsAsFactors=FALSE)

### Convert binary matrix
new=matrix(0, nrow=nrow(my.data), ncol=5)
colnames(new) = colnames(my.data)[10:14]
for (i in 1:nrow(my.data)) {
    r = as.vector(my.data[i,], mode="character")
    g = sapply(strsplit(r[10:length(r)], split=":"), `[[`, 1)
    g[g=="."] = NA
    new[i,] = as.numeric(g)
}
snp=t(new)
stree = nj(dist.gene(snp))
pdf("~/genomics-workshop/phylo.pdf", width=8, height=11)
plot(stree)
dev.off()
