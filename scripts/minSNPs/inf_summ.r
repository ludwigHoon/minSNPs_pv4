library(minSNPs)
library(data.table)

#results <- list.files(pattern = "*.tsv")

f_result_mcc_multi_country <- c("mcc_multi_country.tsv", "mcc_multi_country2.tsv")
f_result_mcc_multi_region <- c("mcc_multi_region.tsv")

f_result_simpson_group_country <- c("simpson_by_group_country.tsv", "simpson_by_group_country2.tsv")
f_result_simpson_group_region <- c("simpson_by_group_region.tsv", "simpson_by_group_region2.tsv")

result_mcc_esea <- process_result_file("mcc_esea.tsv")
result_percent_esea <- process_result_file("percent_esea.tsv")
result_mcc_vietnam <- process_result_file("mcc_vietnam.tsv")
result_percent_vietnam <- process_result_file("percent_vietnam.tsv")

result_simpson <- process_result_file("simpson_old.tsv")


result_m_c <- unlist(lapply(f_result_mcc_multi_country, process_result_file), recursive = FALSE)
result_m_r <- unlist(lapply(f_result_mcc_multi_region, process_result_file), recursive = FALSE)
result_s_c <- unlist(lapply(f_result_simpson_group_country, process_result_file), recursive = FALSE)
result_s_r <- unlist(lapply(f_result_simpson_group_region, process_result_file), recursive = FALSE)


metadata <- fread("metadata.csv")
colnames(metadata)[which(colnames(metadata) == "sample")] <- "isolate"

all_seqs <- read_fasta("simgsk.fasta")
training_seqs <- all_seqs[metadata[partition == "training"]$isolate]
validation_seqs <- all_seqs[metadata[partition == "validation"]$isolate]

###################################### START HERE ######################################
priority <- data.frame(target = c("GOI", "NON_GOI"), priority = c(1, 2))
## ESEA
metadata$is_goi <- ifelse(metadata$region == "ESEA", "GOI", "NON_GOI")
colnames(metadata)[which(colnames(metadata) == "is_goi")] <- "target"
training_metadata <- metadata[partition == "training"]
validation_metadata <- metadata[partition == "validation"]
mcc_esea_result <- summarise_result(result_mcc_esea, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
mcc_esea_result$r_type <- "mcc_esea"
percent_ease_result <- summarise_result(result_percent_esea, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
percent_ease_result$r_type <- "percent_esea"

r_mcc_esea_result <- summarise_result(result_mcc_esea, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE, return_all_intermediate = TRUE)
r_mcc_esea_result$r_type <- "mcc_esea"
r_percent_ease_result <- summarise_result(result_percent_esea, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE, return_all_intermediate = TRUE)
r_percent_ease_result$r_type <- "percent_esea"


mcc_esea_result_sensitive <- summarise_result(result_mcc_esea, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_percent = TRUE)
mcc_esea_result_sensitive$r_type <- "mcc_esea_sensitive"
percent_ease_result_sensitive <- summarise_result(result_percent_esea, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_percent = TRUE)
percent_ease_result_sensitive$r_type <- "percent_esea_sensitive"

r_mcc_esea_result_sensitive <- summarise_result(result_mcc_esea, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE, return_all_intermediate = TRUE, is_percent = TRUE)
r_mcc_esea_result_sensitive$r_type <- "mcc_esea_sensitive"
r_percent_ease_result_sensitive <- summarise_result(result_percent_esea, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE, return_all_intermediate = TRUE, is_percent = TRUE)
r_percent_ease_result_sensitive$r_type <- "percent_esea_sensitive"



metadata[,c("target"):= NULL ]

## Vietnam
metadata$is_goi <- ifelse(metadata$country == "Vietnam", "GOI", "NON_GOI")
colnames(metadata)[which(colnames(metadata) == "is_goi")] <- "target"
training_metadata <- metadata[partition == "training"]
validation_metadata <- metadata[partition == "validation"]

mcc_vietnam_result <- summarise_result(result_mcc_vietnam, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
mcc_vietnam_result$r_type <- "mcc_vietnam"
percent_vietnam_result <- summarise_result(result_percent_vietnam, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
percent_vietnam_result$r_type <- "percent_vietnam"

r_mcc_vietnam_result <- summarise_result(result_mcc_vietnam, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE, return_all_intermediate = TRUE)
r_mcc_vietnam_result$r_type <- "mcc_vietnam"
r_percent_vietnam_result <- summarise_result(result_percent_vietnam, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE, return_all_intermediate = TRUE)
r_percent_vietnam_result$r_type <- "percent_vietnam"


mcc_vietnam_result_sensitive <- summarise_result(result_mcc_vietnam, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_percent = TRUE)
mcc_vietnam_result_sensitive$r_type <- "mcc_vietnam_sensitive"
percent_vietnam_result_sensitive <- summarise_result(result_percent_vietnam, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_percent = TRUE)
percent_vietnam_result_sensitive$r_type <- "percent_vietnam_sensitive"

r_mcc_vietnam_result_sensitive <- summarise_result(result_mcc_vietnam, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE, return_all_intermediate = TRUE, is_percent = TRUE)
r_mcc_vietnam_result_sensitive$r_type <- "mcc_vietnam_sensitive"
r_percent_vietnam_result_sensitive <- summarise_result(result_percent_vietnam, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE, return_all_intermediate = TRUE, is_percent = TRUE)
r_percent_vietnam_result_sensitive$r_type <- "percent_vietnam_sensitive"




all_result <- rbindlist(list(mcc_esea_result, percent_ease_result, mcc_vietnam_result, percent_vietnam_result, mcc_esea_result_sensitive, percent_ease_result_sensitive, mcc_vietnam_result_sensitive, percent_vietnam_result_sensitive))
fwrite(all_result, "mcc_v_percent_esea_vn.csv", row.names = FALSE)
all_r_result <- rbindlist(list(r_mcc_esea_result, r_percent_ease_result, r_mcc_vietnam_result, r_percent_vietnam_result, r_mcc_esea_result_sensitive, r_percent_ease_result_sensitive, r_mcc_vietnam_result_sensitive, r_percent_vietnam_result_sensitive))
fwrite(all_r_result, "mcc_v_percent_esea_vn_raw.csv", row.names = FALSE)


metadata[,c("target"):= NULL ]

###################################### END HERE ######################################



## Region
colnames(metadata)[which(colnames(metadata) == "region")] <- "target"
training_metadata <- metadata[partition == "training"]
validation_metadata <- metadata[partition == "validation"]
priority <- generate_prioritisation(training_metadata[,c("isolate", "target")])

mcc_region_result <- summarise_result(result_m_r, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
mcc_region_result$r_type <- "mcc_multi_region"
simpson_region_result <- summarise_result(result_s_r, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
simpson_region_result$r_type <- "simpson_group_region"
simpson_result_region <- summarise_result(result_simpson, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
simpson_result_region$r_type <- "simpson"

r_mcc_region_result <- summarise_result(result_m_r, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE)
r_mcc_region_result$r_type <- "mcc_multi_region"
r_simpson_region_result <- summarise_result(result_s_r, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE)
r_simpson_region_result$r_type <- "simpson_group_region"
r_simpson_result_region <- summarise_result(result_simpson, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE)
r_simpson_result_region$r_type <- "simpson_region"

colnames(metadata)[which(colnames(metadata) == "target")] <- "region"

## Country
colnames(metadata)[which(colnames(metadata) == "country")] <- "target"
training_metadata <- metadata[partition == "training"]
validation_metadata <- metadata[partition == "validation"]
priority <- generate_prioritisation(training_metadata[,c("isolate", "target")])

mcc_country_result <- summarise_result(result_m_c, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
mcc_country_result$r_type <- "mcc_multi_country"
simpson_country_result <- summarise_result(result_s_c, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
simpson_country_result$r_type <- "simpson_group_country"
simpson_result_country <- summarise_result(result_simpson, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
simpson_result_country$r_type <- "simpson_country"

r_mcc_country_result <- summarise_result(result_m_c, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE)
r_mcc_country_result$r_type <- "mcc_multi_country"
r_simpson_country_result <- summarise_result(result_s_c, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE)
r_simpson_country_result$r_type <- "simpson_group_country"
r_simpson_result_country <- summarise_result(result_simpson, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE)
r_simpson_result_country$r_type <- "simpson_region"

colnames(metadata)[which(colnames(metadata) == "target")] <- "country"

## Can result from country be used to infer region?  i.e., if country is wrong, would the region at least be correct? 
colnames(metadata)[which(colnames(metadata) == "region")] <- "target"
training_metadata <- metadata[partition == "training"]
validation_metadata <- metadata[partition == "validation"]
priority <- generate_prioritisation(training_metadata[,c("isolate", "target")])

mcc_country_2_region_result <- summarise_result(result_m_c, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
mcc_country_2_region_result$r_type <- "mcc_multi_country_2_region"
simpson_country_2_region_result <- summarise_result(result_s_c, training_seqs, validation_seqs, training_metadata, validation_metadata, priority)
simpson_country_2_region_result$r_type <- "simpson_group_country_2_region"


r_mcc_country_2_region_result <- summarise_result(result_m_c, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE)
r_mcc_country_2_region_result$r_type <- "mcc_multi_country_2_region"
r_simpson_country_2_region_result <- summarise_result(result_s_c, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = FALSE)
r_simpson_country_2_region_result$r_type <- "simpson_group_country_2_region"

colnames(metadata)[which(colnames(metadata) == "target")] <- "region"

all_result <- rbindlist(list(mcc_region_result,simpson_region_result, mcc_country_result, simpson_country_result, mcc_country_2_region_result, simpson_country_2_region_result, simpson_result_region, simpson_result_country))
fwrite(all_result, "mcc_simpson_by_group_v_simpson_result.csv", row.names = FALSE)

all_r_region_result <- rbindlist(list(r_mcc_region_result, r_simpson_region_result, r_mcc_country_2_region_result, r_simpson_country_2_region_result, r_simpson_result_region))
all_r_country_result <- rbindlist(list(r_mcc_country_result, r_simpson_country_result, r_simpson_result_country))
fwrite(all_r_region_result, "mcc_simpson_by_group_v_simpson_result_raw_region.csv", row.names = FALSE)
fwrite(all_r_country_result, "mcc_simpson_by_group_v_simpson_result_raw_country.csv", row.names = FALSE)