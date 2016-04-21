require(stringdist)
require(gtools)

medicaments <- read.table('data/medications/Codi_ATC.csv', sep=';', header = TRUE, stringsAsFactors = FALSE)
gcat <- read.csv('data/phenotypes/data.csv', header=TRUE, stringsAsFactors = FALSE)

text <- colnames(gcat)[grep("MEDICACION_.*OTROS", colnames(gcat))]
medi <- colnames(gcat)[grep("MEDICACION_.*_TIPO", colnames(gcat))]

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

