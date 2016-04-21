library(PheWAS)

example <- function() {
  example(PheWAS)
}

exampleUnwind <- function() {

  set.seed(1)
  ex=generateExample()
  
  id.icd9.count=ex$id.icd9.count
  genotypes=ex$genotypes
  phenotypes=createPhewasTable(id.icd9.count)
  results=phewas(phenotypes,genotypes,cores=1, significance.threshold=c("bonferroni"))
  
  phewasManhattan(results, annotate.angle=0, title="My Example PheWAS Manhattan Plot")
}

makeBED <- function() {
  system('plink --noweb --file /imppc/labs/share/lslab/data/Biobank_MEGA_EX/160308_MEGAFinal_Pl_1-13/PLINK_TOP/160308_CF_MEGAFinal_Pl_1-13 --make-bed --out data/GCAT')
}