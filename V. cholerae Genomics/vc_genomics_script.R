#import libraries
library(taxonomizr)
library(rBLAST)
library(filesstrings)

#required files
## Annotation file for genomes -- vcGenomes
## Sequence annotation file from NCBI database #optional


#set working directory
suseWorkingDirectory <- file.choose()

#import metadata
vcGenomes <- read.delim("Annotations_dir.txt")
vcGenomes[is.na(vcGenomes)] <- 0


#qualitative test to check for complete Annotations
sul1Anot <- c()
sul2Anot <- c()
sul3Anot <- c()
dfr1Anot <- c()
dncvAnot <- c()

for (i in 1:length(vcGenomes$NCBI_Strain_ID)) {
  
  if ('sul1' %in% anotFile$Locus) { #test for sul1
    sul1Anot <- c(sul1Anot, 'Yes')
  }else{
    sul1Anot <- c(sul1Anot, 'No')
  }
  if ('sul2' %in% anotFile$Locus) { #test for sul2
    sul2Anot <- c(sul2Anot, 'Yes')
  }else{
    sul2Anot <- c(sul2Anot, 'No')
  }
  if ('sul3' %in% anotFile$Locus) { #test for sul3
    sul3Anot <- c(sul3Anot, 'Yes')
  }else{
    sul3Anot <- c(sul3Anot, 'No')
  }
  if ('dfrA1' %in% anotFile$Locus) { #test for dfrA1
    dfr1Anot <- c(dfr1Anot, 'Yes')
  }else{
    dfr1Anot <- c(dfr1Anot, 'No')
  }
  if ('dncV' %in% anotFile$Locus) { #test for dncV
    dncvAnot <- c(dncvAnot, 'Yes')
  }else{
    dncvAnot <- c(dncvAnot, 'No')
  }
  
}

#append to vcGenomes
vcGenomes$sul1ANOT <- sul1Anot
vcGenomes$sul2ANOT <- sul2Anot
vcGenomes$sul3ANOT <- sul3Anot
vcGenomes$dfrANOT <- dfr1Anot
vcGenomes$dncVANOT <- dncvAnot


#BLAST
#move files to individual files to create blast dbs
#do this once
for (fileN in list.files("Genome_Prep/")) {
  myFiles <- strsplit(fileN, ".fna")[[1]][1]
  dir.create(paste0("FullGenomes/", myFiles))
  file.move(paste0("Genome_Prep/", fileN), paste0("FullGenomes/", myFiles))
}

#implement blast
blastResult <- list()

#input queries
sul1Dna <- readDNAStringSet("Queries/sul1.fasta")
sul2DNA <- readDNAStringSet("Queries/sul2.fasta")
sul3DNA <- readDNAStringSet("Queries/sul3.fasta") #only available for E.coli
dfrA1 <- readDNAStringSet("Queries/dfrA1.fasta") #renamed from dhfrII
dncvDNA <- readDNAStringSet("Queries/dncV.fasta")  #GOTTEN FROM https://www.genome.jp/entry/-f+-n+n+vch:VC_0179
cap2DNA <- readDNAStringSet("Queries/cap2.fasta")
cap3DNA <- readDNAStringSet("Queries/cap3.fasta")
capVDNA <- readDNAStringSet("Queries/capV.fasta")
vsp1DNA <- readDNAStringSet("Queries/vsp1.fasta")

#accept percent.ident
sul1Result <- c()
sul2Result <- c()
sul3Result <- c()
dfra1Result <- c()
dncvResult <- c()
cap2Result <- c()
cap3Result <- c()
capVResult <- c()
strainId <- c()
vsp1Result <- c()

#setup reusable script for BLAST-ing
for (iGen in vcGenomes$NCBI_Strain_ID) {
  print(which(iGen == vcGenomes$NCBI_Strain_ID))
  
  strainId <- c(strainId, iGen)
  
  #implement blast
  makeblastdb(paste0(suseWorkingDirectory, "/FullGenomes/", iGen, paste0("/",iGen, '.fna')))
  currBl <- blast(db = paste0(suseWorkingDirectory, "/FullGenomes/", iGen, paste0("/",iGen, '.fna')))
  
  #sul1Cl
  sul1Cl <- predict(currBl, sul1Dna)
  sul1Result <- c(sul1Result, list(sul1Cl))
  
  #sul2Cl
  sul2Cl <- predict(currBl, sul2DNA)
  sul2Result <- c(sul2Result, list(sul2Cl))
  
  #sul3Cl
  sul3Cl <- predict(currBl, sul3DNA)
  sul3Result <- c(sul3Result, list(sul3Cl))
  
  #dfra1Cl
  dfraCl <- predict(currBl, dfrA1)
  dfra1Result <- c(dfra1Result, list(dfraCl))
  
  #dncvCl
  dncvCl <- predict(currBl, dncvDNA)
  dncvResult <- c(dncvResult, list(dncvCl))
  
  #cap2Cl
  cap2Cl <- predict(currBl, cap2DNA)
  cap2Result <- c(cap2Result, list(cap2Cl))
  
  #cap3Cl
  cap3Cl <- predict(currBl, cap3DNA)
  cap3Result <- c(cap3Result, list(cap3Cl))
  
  #capVCl
  capVCl <- predict(currBl, capVDNA)
  capVResult <- c(capVResult, list(capVCl))
  
  #vsp1CL
  vsp1Cl <- predict(currBl, vsp1DNA)
  vsp1Result <- c(vsp1Result, list(vsp1Cl))
  
}

#process all *Results

sul1ResultProc <- c()
sul2ResultProc <- c()
sul3ResultProc <- c()
dncvResultProc <- c()
cap2ResultProc <- c()
cap3ResultProc <- c()
capVResultProc <- c()
vsp1ResultProc <- c()

for (i in 1:length(sul1Result)) {
  if(class(sul1Result[[i]]$Perc.Ident) == 'logical'){
    sul1ResultProc <- c(sul1ResultProc, 0)
  }else{
    sul1ResultProc <- c(sul1ResultProc, sul1Result[[i]]$Perc.Ident[1])
  }
  
  if(class(sul2Result[[i]]$Perc.Ident) == 'logical'){
    sul2ResultProc <- c(sul2ResultProc, 0)
  }else{
    sul2ResultProc <- c(sul2ResultProc, sul2Result[[i]]$Perc.Ident[1])
  }
  
  if(class(sul3Result[[i]]$Perc.Ident) == 'logical'){
    sul3ResultProc <- c(sul3ResultProc, 0)
  }else{
    sul3ResultProc <- c(sul3ResultProc, sul3Result[[i]]$Perc.Ident[1])
  }
  
  if(class(dncvResult[[i]]$Perc.Ident) == 'logical'){
    dncvResultProc <- c(dncvResultProc, 0)
  }else{
    dncvResultProc <- c(dncvResultProc, dncvResult[[i]]$Perc.Ident[1])
  }
  
  if(class(cap2Result[[i]]$Perc.Ident) == 'logical'){
    cap2ResultProc <- c(cap2ResultProc, 0)
  }else{
    cap2ResultProc <- c(cap2ResultProc, cap2Result[[i]]$Perc.Ident[1])
  }
  
  if(class(cap3Result[[i]]$Perc.Ident) == 'logical'){
    cap3ResultProc <- c(cap3ResultProc, 0)
  }else{
    cap3ResultProc <- c(cap3ResultProc, cap3Result[[i]]$Perc.Ident[1])
  }
  
  if(class(capVResult[[i]]$Perc.Ident) == 'logical'){
    capVResultProc <- c(capVResultProc, 0)
  }else{
    capVResultProc <- c(capVResultProc, capVResult[[i]]$Perc.Ident[1])
  }
  
  if(class(vsp1Result[[i]]$Perc.Ident) == 'logical'){
    vsp1ResultProc <- c(vsp1ResultProc, 0)
  }else{
    vsp1ResultProc <- c(vsp1ResultProc, vsp1Result[[i]]$Perc.Ident[1])
  }
}

length(sul1ResultProc)
length(sul2ResultProc)
length(sul3ResultProc)
length(dncvResultProc)
length(cap2ResultProc)
length(cap3ResultProc)
length(capVResultProc)
length(vsp1ResultProc)

vcGenomes$sul1BLAST <- sul1ResultProc
vcGenomes$sul2BLAST <- sul2ResultProc
vcGenomes$sul3BLAST <- sul3ResultProc
#vcGenomes$dfraBLAST <- 
vcGenomes$dncVBLAST <- dncvResultProc



#FINAL PROCESSED DATA
finalProcessedData <- data.frame("BIOCYC_NAME" = vcGenomes$BIOCYC_DB_NAME,
                                 "BIOCYC_DB_NAME" = strainId,
                                 "sul1BLAST" = vcGenomes$sul1BLAST,
                                 'sul2BLAST' = vcGenomes$sul2BLAST,
                                 "sul3BLAST" = vcGenomes$sul3BLAST,
                                 "dncVBLAST" = vcGenomes$dncVBLAST,
                                 "cap2" = cap2ResultProc,
                                 "cap3" = cap3ResultProc,
                                 "capV" = capVResultProc,
                                 "vsp1" = vsp1ResultProc)

View(finalProcessedData)


#collecting counts
cbassResStatus <- c()

for (bug in 1:nrow(finalProcessedData)) {
  currBug <- finalProcessedData[bug,]
  
  if ((currBug$sul1BLAST == 100 | currBug$sul2BLAST == 100) & (currBug$dncVBLAST == 100 & currBug$cap2 > 90 & currBug$cap3 == 100 & currBug$capV == 100) ) {
    cbassResStatus <- c(cbassResStatus, '100_100')
    
  }else if ((currBug$sul1BLAST == 100 | currBug$sul2BLAST >= 99) & (currBug$dncVBLAST==0 & currBug$cap2==0 & currBug$cap3==0 & currBug$capV==0) ) {
    cbassResStatus <- c(cbassResStatus, '0_100')
    #print(currBug)
  }else if ((currBug$sul1BLAST == 0 & currBug$sul2BLAST == 0) & (currBug$dncVBLAST >= 99 & currBug$cap2 > 90 & currBug$cap3 == 100 & currBug$capV == 100) ) {
    cbassResStatus <- c(cbassResStatus, '100_0')
    #print(currBug)
  }else if ((currBug$sul1BLAST == 0 & currBug$sul2BLAST == 0) & (currBug$dncVBLAST==0 & currBug$cap2==0 & currBug$cap3==0 & currBug$capV==0) ) {
    cbassResStatus <- c(cbassResStatus, '0_0')
  }else{
    print(currBug)
  }
  
}

#plotting
table(cbassResStatus)

barplot(table(cbassResStatus), col = myColorList$Orange)
legend("topleft", legend = paste0("N =",length(cbassResStatus)), col = "#F6512F", pch = 19)

#pandemic status
pandemicGen <- c()
for (yearN in vcGenomes$Collection_Year) {
  if (yearN == 0) {
    pandemicGen <- c(pandemicGen, "Unknown")
  }else if (yearN >= 1961 ) {
    #print(yearN)
    pandemicGen <- c(pandemicGen, "7th Gen")
  }else if(yearN >= 1899 & yearN <= 1923){
    #print(yearN)
    pandemicGen <- c(pandemicGen, "6th Gen")
  }else{
    print(yearN)
    pandemicGen <- c(pandemicGen, "pre-6th Gen")
  }
}


#append pandemic and resistance status to the finalProcessedData
finalProcessedData$cbassResStatus <- cbassResStatus
finalProcessedData$pandemicStatus <- pandemicGen


#processing cbass status in each genome -- qualitative

CBASS_PRESENT <- c()
SUL_GENE_PRESENT <- c()
VSP1_LOCUS_PRESENT <- c()

for (i in 1:nrow(finalProcessedData)) {
  if (finalProcessedData[i,]$cbassResStatus == "100_0" | finalProcessedData[i,]$cbassResStatus == "100_100") {
    CBASS_PRESENT <- c(CBASS_PRESENT, "Yes") 
  }else{
    CBASS_PRESENT <- c(CBASS_PRESENT, "No")
  }
  
  if (finalProcessedData[i,]$cbassResStatus == "0_100" | finalProcessedData[i,]$cbassResStatus == "100_100") {
    SUL_GENE_PRESENT <- c(SUL_GENE_PRESENT, "Yes") 
  }else{
    SUL_GENE_PRESENT <- c(SUL_GENE_PRESENT, "No")
  }
  
  if (finalProcessedData[i,]$vsp1 == 0.00) {
    VSP1_LOCUS_PRESENT <- c(VSP1_LOCUS_PRESENT, "No") 
  }else{
    VSP1_LOCUS_PRESENT <- c(VSP1_LOCUS_PRESENT, "Yes")
  }
  
}


CBASS_PRESENT_HT <- c()
SUL_GENE_PRESENT_HT <- c()
VSP1_LOCUS_PRESENT_HT <- c()

for (i in 1:nrow(finalProcessedData)) {
  if (finalProcessedData[i,]$cbassResStatus == "100_0" | finalProcessedData[i,]$cbassResStatus == "100_100") {
    CBASS_PRESENT_HT <- c(CBASS_PRESENT_HT, 1) 
  }else{
    CBASS_PRESENT_HT <- c(CBASS_PRESENT_HT, 0)
  }
  
  if (finalProcessedData[i,]$cbassResStatus == "0_100" | finalProcessedData[i,]$cbassResStatus == "100_100") {
    SUL_GENE_PRESENT_HT <- c(SUL_GENE_PRESENT_HT, 1) 
  }else{
    SUL_GENE_PRESENT_HT <- c(SUL_GENE_PRESENT_HT, 0)
  }
  
  if (finalProcessedData[i,]$vsp1 == 0.00) {
    VSP1_LOCUS_PRESENT_HT <- c(VSP1_LOCUS_PRESENT_HT, 0) 
  }else{
    VSP1_LOCUS_PRESENT_HT <- c(VSP1_LOCUS_PRESENT_HT, 1)
  }
  
}

#append cbass status to the finalProcessedData
finalProcessedData$CBASS_PRESENT <- CBASS_PRESENT
finalProcessedData$SUL_PRESENT <- SUL_GENE_PRESENT
finalProcessedData$VSP1_PRESENT <- VSP1_LOCUS_PRESENT


#Fisher testing
table(finalProcessedData$cbassResStatus)
fis.tab <- data.frame("cbassPos" = c(43, 15), "cbassNeg" = c(15, 23))
row.names(fis.tab) <- c("sulPos", "sulNeg")

fis.tab
cbass.fish <- fisher.test(fis.tab)
cbass.fish$p.value

#plotting
pdf("report.pdf", width = 80, height = 120, paper = 'a4')
#bar chart for abundance of cbass operon
par(mfrow = c(2,2), mar = c(2,2,2,2))
barplot(table(cbassResStatus), col = myColorList$Orange, ylim = c(0,50))
legend("topleft", legend = paste0("N =",length(cbassResStatus)), col = "#F6512F", pch = 19)

#pie chart for individuals genotypes and the abundance of the pandemics
for (i in unique(finalProcessedData$cbassResStatus)) {
  cb <- subset(finalProcessedData, cbassResStatus == i)
  pie(table(cb$pandemicStatus), main = i, col =  c(myColorList$greens[5], myColorList$blues[5], myColorList$Red[5]))
}

#pie chart for individuals genotypes and the abundance of VSP
for (i in unique(finalProcessedData$cbassResStatus)) {
  cb <- subset(finalProcessedData, cbassResStatus == i)
  pie(table(cb$VSP1_PRESENT), main = i, col =  c(myColorList$greens[5], myColorList$blues[5], myColorList$Red[5]))
}

#mosiac plot of cbass and sul-resistance gene counts across all genomes
mosaicplot(fis.tab,
           main = 'Contingency', 
           shade = F, 
           color = myColorList$greens[5], 
           border = myColorList$greens[5], 
           las = 1)

dev.off()

