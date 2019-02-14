library(stringr)
library(stringi)
library(readxl)

MotifPhos <- read_excel("5198_Phos_for_MotifX.xlsx", 1)
Fasta_file <- read_excel("5198_Fasta.xlsx", 1)
file_prefix <- "5198_Phos_020618_Motif_Final.xlsx"

phosmotif<- "[KR][KR].[ST]\\*."

MotifPhos$Accession <- gsub(";.*", "", MotifPhos$Accession)
MotifPhos <- subset(MotifPhos, Accession != 1000 )
MotifPhos$PhosOnly2 <- str_extract(MotifPhos$Modifications, "P.+\\]")
MotifPhos$PhosOnly3 <- str_extract_all(MotifPhos$PhosOnly2, "[STY]\\d+")
MotifPhos$PhosOnly3[MotifPhos$PhosOnly3=="character(0)"] <- NA
MotifPhos$PhosOnly4 <- str_extract_all(MotifPhos$PhosOnly3, "\\d+")
MotifPhos$PhosOnly4[MotifPhos$PhosOnly4=="character(0)"] <- 0
MotifPhos$PhosOnly6 <- str_extract_all(MotifPhos$PhosOnly3, "[STY]")
MotifPhos$PhosOnly6[MotifPhos$PhosOnly6=="character(0)"] <- NA
MotifPhos$Seq2 <- MotifPhos$Sequence

    
Fasta_file <-subset(Fasta_file, Fasta_file$Accession %in% MotifPhos$Accession)
test_merge <- merge(MotifPhos, Fasta_file, "Accession", all.x = FALSE)
test_merge$startAA <- str_locate(test_merge$Sequence.y, test_merge$Sequence.x)
test_merge$endAA <- test_merge$startAA + nchar(test_merge$Sequence.x)
test_merge$addStart <- substr(test_merge$Sequence.y, test_merge$startAA-7, test_merge$startAA-1)
test_merge$addEnd <- substr(test_merge$Sequence.y, test_merge$endAA+1, test_merge$endAA+7)
test_merge$ext_AA <- ""


for (i in 1:nrow(test_merge)) {
  phos_res <- rev(unlist(test_merge$PhosOnly4[i]))
  aa_list <- unlist(test_merge$PhosOnly6[i])
  loc_list <- as.numeric(phos_res) + test_merge$startAA[i]-1
  df <- data.frame(cbind(aa_list, rev(loc_list)))
  colnames(df)<-c("aa", "loc")
  df$combo <- str_c(df$aa, df$loc)
  test_merge$ext_AA[i] <- list(df$combo)
  for (j in 1:length(phos_res)) {
    set_loc <- as.numeric(phos_res[j])
    if (set_loc > 0) {stri_sub(test_merge$Seq2[i], set_loc+1, set_loc) <- "*"} 
  }
}


test_merge$extended <- str_c(test_merge$addStart,test_merge$Seq2, test_merge$addEnd)

test_merge$testmotif <- ""
test_merge$testmotif <- grepl(phosmotif, test_merge$extended)


final_phos <- data.frame(cbind(test_merge$Accession, test_merge$Description, test_merge$Sequence.x, test_merge$Modifications, 
                    test_merge$Seq2, test_merge$ext_AA, test_merge$extended, test_merge$testmotif))


col_headers <- c("Accession", "Description", "Sequence", "Modifications", "MotifX", "ProteinAA", "Extended Motif", "Test Motif")
colnames(final_phos) <- col_headers


wb <- createWorkbook()
addWorksheet(wb, deparse(substitute(final_phos)))
writeData(wb, sheet =1, final_phos)  
saveWorkbook(wb, file_prefix, overwrite = TRUE)

