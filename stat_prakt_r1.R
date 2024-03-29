### Einlesen Dateien
FIRE_3_clin_data <- read_excel("Desktop/Statistisches Praktikum/FIRE-3_clin_data.xlsx")
FIRE_3_EMT_genes_final <- read_excel("Desktop/Statistisches Praktikum/FIRE-3_EMT_genes_final.xlsx")

### Umbenennen von CRF in pat_nr fuer merge
FIRE_3_EMT_genes_final$pat_nr <- FIRE_3_EMT_genes_final$CRF   

str(FIRE_3_clin_data)
str(FIRE_3_EMT_genes_final)

### merge by pat_nr, NA aussortiert
fire3_complete <- merge(FIRE_3_clin_data, FIRE_3_EMT_genes_final, by = "pat_nr")

str(FIRE_3_EMT_genes_final)
