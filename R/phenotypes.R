library(PheWAS)
library(stringdist)
library(XLConnect)
library(gtools)

getParticipants <- function() {
  family <- read.table('data/genotypes/GCAT.fam', header = FALSE)
  participants <- read.table('data/phenotypes/participants.csv')
  participants <- participants[participants$Sample.name %in% family$V2,]
  participants <- participants[, c('Sample.Id', 'Sample.name')]
  participants
}

getDiseases <- function(participants) {

  gcat <- read.csv('data/phenotypes/data.csv', header=TRUE, stringsAsFactors = FALSE)
  phenotypes_description <- read.csv('data/phen_description.csv', header=TRUE, stringsAsFactors = FALSE)
  phenotypes <- gcat[, append(c('entity_id'), phenotypes_description$phenotype)]
  phenotypes$entity_id <- gsub('=', '', phenotypes$entity_id)
  phenotypes <- cbind(phenotypes[1], apply(phenotypes[2:ncol(phenotypes)], 2, function(column) {
    ifelse(column == 'false', F, T)
  }))
  
  write.csv(phenotypes, 'output/diseases.csv', row.names = FALSE)

  names(phenotypes) <- append(c('entity_id'), phenotypes_description$icd9)
  phenotypes <- merge(phenotypes, participants, by.x = 'entity_id', by.y = 'Sample.Id')
  phenotypes <- phenotypes[, !(colnames(phenotypes) %in% c('entity_id'))]
  names(phenotypes)[which(colnames(phenotypes) == 'Sample.name')] <- 'id'
  phenotypes
}

getGenotypes <- function(outcome, threshold, types = NA) {
  associations <- read.csv(sprintf('data/gwas-association-%s.tsv', outcome), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  associations <- associations[which(associations$PVALUE_MLOG > threshold),]
  if (!is.na(types)) {
    associations <- subset(associations, associations$CONTEXT %in% types)
  }
  write.table(associations$SNPS, 'output/snps.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
  system('plink --noweb --recodeA --bfile data/genotypes/GCAT --extract output/snps.txt --out output/r_genotypes')
  genotypes <- read.table("output/r_genotypes.raw", header=TRUE)[,c(-1,-3:-6)]
  names(genotypes)[1] = "id"
  genotypes
}

do <- function() {

  genotypes <- getGenotypes('hypertension', 10, c('missense_variants'))
  participants <- getParticipants()
  diseases <- getDiseases(participants) 
  
  result <- phewas(phenotypes = diseases, genotypes = genotypes, cores=1, significance.threshold=c("bonferroni"))
  phewasManhattan(result, annotate.angle=0, title="My Example PheWAS Manhattan Plot", annotate.phenotype = TRUE, annotate.snp = TRUE)
}