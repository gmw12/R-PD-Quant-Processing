MotifPhos <- read_excel("5081_GluC_Ab_Phos_for_Motif.xlsx", 1)
Fasta_file <- read_excel("5081_Fasta.xlsx", 1)
file_prefix <- "5081_GluC_Ab_Motif_Final.xlsx"

phosmotif<- "[KR][KR].[ST]\\*."

MotifPhos$Accession <- gsub(";.*", "", MotifPhos$Accession)
MotifPhos <- subset(MotifPhos, Accession != 1000 )


MotifPhos$PhosOnly <- substring(MotifPhos$Modifications, regexpr("Phospho", MotifPhos$Modifications))
MotifPhos$PhosOnly <- substring(MotifPhos$PhosOnly, regexpr("\\[", MotifPhos$PhosOnly))
MotifPhos$PhosOnly <- str_replace(MotifPhos$PhosOnly, "\\[", "")
MotifPhos$PhosOnly <- str_replace(MotifPhos$PhosOnly, "\\]", "")
MotifPhos$PhosOnly <- gsub("\\(.*?\\)", "", MotifPhos$PhosOnly)

MotifPhos$AA1 <- MotifPhos$PhosOnly
MotifPhos$AA1 <- str_replace_all(MotifPhos$AA1," ", "")
MotifPhos$AA1 <- str_replace_all(MotifPhos$AA1,"././.", "")
MotifPhos$AA1 <- str_replace_all(MotifPhos$AA1,"./.", "")
MotifPhos$AA1 <- gsub("\\d*", "", MotifPhos$AA1)

MotifPhos$AA2 <- gsub("^.", "", MotifPhos$AA1)
MotifPhos$AA2 <- gsub("^;", "", MotifPhos$AA2)

MotifPhos$AA3 <- gsub("^.", "", MotifPhos$AA2)
MotifPhos$AA3 <- gsub("^;", "", MotifPhos$AA3)

MotifPhos$AA4 <- gsub("^.", "", MotifPhos$AA3)
MotifPhos$AA4 <- gsub("^;", "", MotifPhos$AA4)

MotifPhos$AA5 <- gsub("^.", "", MotifPhos$AA4)
MotifPhos$AA5 <- gsub("^;", "", MotifPhos$AA5)


MotifPhos$AA1 <- gsub(";.*", "", MotifPhos$AA1)
MotifPhos$AA2 <- gsub(";.*", "", MotifPhos$AA2)
MotifPhos$AA3 <- gsub(";.*", "", MotifPhos$AA3)
MotifPhos$AA4 <- gsub(";.*", "", MotifPhos$AA4)
MotifPhos$AA5 <- gsub(";.*", "", MotifPhos$AA5)

MotifPhos$PhosMark1 <- str_replace_all(MotifPhos$PhosOnly, "S", "*")
MotifPhos$PhosMark1 <- str_replace_all(MotifPhos$PhosMark1, "T", "*")
MotifPhos$PhosMark1 <- str_replace_all(MotifPhos$PhosMark1, "Y", "*")

MotifPhos$PhosMark1 <- str_replace_all(MotifPhos$PhosMark1, "\\*", "")
MotifPhos$PhosMark1 <- str_replace_all(MotifPhos$PhosMark1, ";", ",")
MotifPhos$PhosMark1 <- str_replace_all(MotifPhos$PhosMark1, "\\/", "")
MotifPhos$PhosMark1 <- gsub(",.$", "", MotifPhos$PhosMark1)

MotifPhos$PhosMark2 <- gsub("^\\d*", "", MotifPhos$PhosMark1)
MotifPhos$PhosMark2 <- gsub("^,", "", MotifPhos$PhosMark2)
MotifPhos$PhosMark2 <- gsub("^ ", "", MotifPhos$PhosMark2)

MotifPhos$PhosMark3 <- gsub("^\\d*", "", MotifPhos$PhosMark2)
MotifPhos$PhosMark3 <- gsub("^,", "", MotifPhos$PhosMark3)
MotifPhos$PhosMark3 <- gsub("^ ", "", MotifPhos$PhosMark3)

MotifPhos$PhosMark4 <- gsub("^\\d*", "", MotifPhos$PhosMark3)
MotifPhos$PhosMark4 <- gsub("^,", "", MotifPhos$PhosMark4)
MotifPhos$PhosMark4 <- gsub("^ ", "", MotifPhos$PhosMark4)

MotifPhos$PhosMark5 <- gsub("^\\d*", "", MotifPhos$PhosMark4)
MotifPhos$PhosMark5 <- gsub("^,", "", MotifPhos$PhosMark5)
MotifPhos$PhosMark5 <- gsub("^ ", "", MotifPhos$PhosMark5)

MotifPhos$PhosMark1 <- gsub(",.*", "", MotifPhos$PhosMark1)
MotifPhos$PhosMark2 <- gsub(",.*", "", MotifPhos$PhosMark2)
MotifPhos$PhosMark3 <- gsub(",.*", "", MotifPhos$PhosMark3)
MotifPhos$PhosMark4 <- gsub(",.*", "", MotifPhos$PhosMark4)
MotifPhos$PhosMark5 <- gsub(",.*", "", MotifPhos$PhosMark5)

MotifPhos$PhosMark1 <- as.numeric(MotifPhos$PhosMark1)
MotifPhos$PhosMark2 <- as.numeric(MotifPhos$PhosMark2)
MotifPhos$PhosMark3 <- as.numeric(MotifPhos$PhosMark3)
MotifPhos$PhosMark4 <- as.numeric(MotifPhos$PhosMark4)
MotifPhos$PhosMark5 <- as.numeric(MotifPhos$PhosMark5)

MotifPhos$PhosMark2[is.na(MotifPhos$PhosMark2)] <- 0
MotifPhos$PhosMark3[is.na(MotifPhos$PhosMark3)] <- 0
MotifPhos$PhosMark4[is.na(MotifPhos$PhosMark4)] <- 0
MotifPhos$PhosMark5[is.na(MotifPhos$PhosMark5)] <- 0


MotifPhos <- MotifPhos[complete.cases(MotifPhos),]

MotifPhos$MotifX <- str_c(substring(MotifPhos$Sequence, 1, MotifPhos$PhosMark1),"*",
                          substring(MotifPhos$Sequence, MotifPhos$PhosMark1+1, nchar(MotifPhos$Sequence)) )



create_motifx <- function(df,x){
  for (i in 1:nrow(MotifPhos))
  {
    if (df[i] > 0)
    {
      MotifPhos$MotifX[i] <<- str_c(substring(MotifPhos$MotifX[i], 1, df[i]+x),"*",
                                    substring(MotifPhos$MotifX[i], df[i]+x+1, nchar(MotifPhos$MotifX[i])))
    }else
    {
      MotifPhos$MotifX[i] <<- MotifPhos$MotifX[i]
    }
  }
}

create_motifx(MotifPhos$PhosMark2, 1)
create_motifx(MotifPhos$PhosMark3, 2)
create_motifx(MotifPhos$PhosMark4, 3)
create_motifx(MotifPhos$PhosMark5, 4)




Fasta_file <-subset(Fasta_file, Fasta_file$Accession %in% MotifPhos$Accession)


test_merge <- merge(MotifPhos, Fasta_file, "Accession", all.x = FALSE)

test_merge$startAA <- str_locate(test_merge$Sequence.y, test_merge$Sequence.x)
test_merge$endAA <- test_merge$startAA + nchar(test_merge$Sequence.x)

test_merge$ProteinMark1<-""
test_merge$ProteinMark2<-""
test_merge$ProteinMark3<-""
test_merge$ProteinMark4<-""
test_merge$ProteinMark5<-""

for (i in 1:nrow(test_merge))
{
  if (test_merge$PhosMark1[i]>0) {test_merge$ProteinMark1[i] <- test_merge$PhosMark1[i] + test_merge$startAA[i] -1}
  if (test_merge$PhosMark2[i]>0) {test_merge$ProteinMark2[i] <- test_merge$PhosMark2[i] + test_merge$startAA[i] -1}  
  if (test_merge$PhosMark3[i]>0) {test_merge$ProteinMark3[i] <- test_merge$PhosMark3[i] + test_merge$startAA[i] -1}
  if (test_merge$PhosMark4[i]>0) {test_merge$ProteinMark4[i] <- test_merge$PhosMark4[i] + test_merge$startAA[i] -1}  
  if (test_merge$PhosMark5[i]>0) {test_merge$ProteinMark5[i] <- test_merge$PhosMark5[i] + test_merge$startAA[i] -1}  
  
}


test_merge$ProteinMark1 <- str_c(test_merge$AA1,  test_merge$ProteinMark1)
test_merge$ProteinMark2 <- str_c(test_merge$AA2,  test_merge$ProteinMark2)
test_merge$ProteinMark3 <- str_c(test_merge$AA3,  test_merge$ProteinMark3)
test_merge$ProteinMark4 <- str_c(test_merge$AA4,  test_merge$ProteinMark4)
test_merge$ProteinMark5 <- str_c(test_merge$AA5,  test_merge$ProteinMark5)

test_merge$ProteinPhosID <- str_c(test_merge$ProteinMark1, ", ", test_merge$ProteinMark2, ", ", 
                                 test_merge$ProteinMark3, ", ", test_merge$ProteinMark4, ", ",
                                 test_merge$ProteinMark5)

test_merge$ProteinPhosID <- gsub("(, )*$", "", test_merge$ProteinPhosID)

test_merge$NtermAA <- ""
test_merge$CtermAA <- ""
test_merge$MotifSearch <-""
for (i in 1:nrow(test_merge))
{
test_merge$NtermAA[i] <- substring(test_merge$Sequence.y[i], test_merge$startAA[i]-5, test_merge$startAA[i]-1)
test_merge$CtermAA[i] <- substring(test_merge$Sequence.y[i], test_merge$endAA[i], test_merge$endAA[i]+4)
test_merge$MotifSearch[i] <- str_c(test_merge$NtermAA[i], test_merge$MotifX[i] ,test_merge$CtermAA[i])
}

test_merge$testmotif <- grepl(phosmotif, test_merge$MotifSearch)


final_phos <- cbind.data.frame(test_merge$Accession, test_merge$Description, test_merge$Sequence.x, test_merge$Modifications, 
                    test_merge$MotifX, test_merge$ProteinPhosID, test_merge$MotifSearch, test_merge$testmotif)

col_headers <- c("Accession", "Description", "Sequence", "Modifications", "MotifX", "ProteinAA", "Extended Motif", "Test Motif")
colnames(final_phos) <- col_headers



wb <- createWorkbook()
addWorksheet(wb, deparse(substitute(final_phos)))
writeData(wb, sheet =1, final_phos)  
saveWorkbook(wb, file_prefix, overwrite = TRUE)
