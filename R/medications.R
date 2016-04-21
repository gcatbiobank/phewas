library(stringdist)
library(gtools)
library(PheWAS)

source('R/phenotypes.R')

getParticipantsATC <- function(participants, granularity = 3) {
  
  medicaments <- read.table('data/medications/Codi_ATC.csv', sep=';', header = TRUE, stringsAsFactors = FALSE)
  gcat <- read.csv('data/phenotypes/data.csv', header=TRUE, stringsAsFactors = FALSE)
  
  text <- colnames(gcat)[grep("MEDICAMENTO_.*OTROS", colnames(gcat))]
  medi <- colnames(gcat)[grep("^MEDICAMENTO_.*", colnames(gcat))]
  
  text <- names(gcat)[grepl("_MEDICACION_.*OTROS", names(gcat))]
  medi <- names(gcat)[grepl('_MEDICACION_.*_TIPO', names(gcat))]
  
  
  gcat.medi <- gcat[, append('entity_id', medi)]
  gcat.medi.reshape <- reshape(gcat.medi, varying = medi, v.names = "MEDICATION", timevar = "CONDITION", times = medi, direction = "long")
  gcat.medi.reshape <- subset(gcat.medi.reshape, gcat.medi.reshape$MEDICATION != '')
  gcat.medi.reshape <- unique(gcat.medi.reshape[, c('entity_id', 'MEDICATION')])
  gcat.medi.reshape$MEDICATION <- gsub('\\.', ';', gcat.medi.reshape$MEDICATION)
  gcat.medi.reshape$MEDICATION <- gsub('\\,', ';', gcat.medi.reshape$MEDICATION)
  
  df <- data.frame(do.call(rbind, apply(gcat.medi.reshape, 1, function(med) {
    sp <- strsplit(med['MEDICATION'], ';')
    c(med['entity_id'], substring(sp[[1]], 1, granularity))
  })), row.names = NULL)
  
  df <- unique(df[, c('entity_id', 'V2')])
  df$entity_id <- gsub('=', '', df$entity_id)
  df$ATC <- T
  
  phenotypes <- reshape(df, timevar='V2', idvar = c('entity_id'), direction = 'wide')
  names(phenotypes) <- gsub('^ATC.', '', names(phenotypes))
  phenotypes[is.na(phenotypes)] <- FALSE
  
  write.csv(phenotypes, 'output/drugs.csv', row.names = FALSE)
  
  phenotypes <- merge(phenotypes, participants, by.x = 'entity_id', by.y = 'Sample.Id')
  phenotypes <- phenotypes[, !(colnames(phenotypes) %in% c('entity_id'))]
  names(phenotypes)[which(colnames(phenotypes) == 'Sample.name')] <- 'id'
  phenotypes$id <- as.character(phenotypes$id)
  
  phenotypes
}

curateParticipantsATC <- function() {
  
  medicaments.text <- do.call(smartbind, lapply(text, function(columna) {
    valors <- unique(gcat[, columna])
    valors <- valors[!is.na(valors) & valors != '']
    df <- data.frame()
    if (length(valors) > 0) {
      df <- data.frame(columna=columna, medicament=valors)
    }
    df
  }))
  
  ### Compute distances
  
  x <- toupper(medicaments.text$medicament)
  table <- as.character(medicaments$NOMBRE)
  
  #' Corregir els medicaments mal escrits
  #' @param x
  #' @param table
  #' @param method
  #' @param maxDist
  correct <- function(x, table, method="osa", maxDist = 2) {
    
    x <- gsub('[[:digit:]]+', '', x)
    x <- gsub('( )*MG', '', x)
    x <- gsub('[[:punct:]]', '', x)
    x <- gsub('[[:space:]]', '', x)
    
    matches <- amatch(x, table, method=method, maxDist=maxDist)
    medicaments.text$CODI_ATC <- NA
    medicaments.text$NOMBRE <- NA
    
    for (item in 1:length(matches)) {
      if (!is.na(matches[item])) {
        medicaments.text[item,]$CODI_ATC <- medicaments[matches[item],]$CODI_ATC
        medicaments.text[item,]$NOMBRE <- medicaments[matches[item],]$NOMBRE
      }
    }
    medicaments.text[, c('CODI_ATC', 'NOMBRE')]
  }
  
  a <- medicaments.text
  a <- cbind(a, correct(x, table, maxDist = 1))
  a <- cbind(a, correct(x, table, maxDist = 2))
  a <- cbind(a, correct(x, table, maxDist = 3))
  a <- cbind(a, correct(x, table, maxDist = 4))
  
  length(which(!is.na(correct(x, table, maxDist = 1)$CODI_ATC)))
  length(which(!is.na(correct(x, table, maxDist = 2)$CODI_ATC)))
  length(which(!is.na(correct(x, table, maxDist = 3)$CODI_ATC)))
}

do <- function() {
  
  genotypes <- getGenotypes('metabolic disease', 20)
  participants <- getParticipants()
  drugs <- getParticipantsATC(participants, 5)
  
  genotypes <- getGenotypes('t1dm', 10, c('missense_variant'))
  genotypes <- getGenotypes('t2dm', 10, c('missense_variant'))
  genotypes <- getGenotypes('cholesterol', 10, c('missense_variant'))
  
  participants <- getParticipants()
  drugs <- getParticipantsATC(participants, 3)
  
  result <- phewas(phenotypes = drugs, genotypes = genotypes, cores=4, significance.threshold=c("bonferroni"))
  result <- subset(result, !is.na(result$p))
  result$value <- result$phenotype
  
  phewasManhattan(result, annotate.angle=0, title="My Example PheWAS Manhattan Plot", annotate.phenotype = TRUE, annotate.snp = TRUE)
}