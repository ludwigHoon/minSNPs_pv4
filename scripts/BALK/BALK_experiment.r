library(minSNPs)
library(data.table)

## Temporary patch for unforeseen issue with using BALK and the SNPs in the data:
## TL;DR: Some of the randomly sampled SNPs have only 1 state, which causes BALK to fail.
## The reason for single state SNP is that the data used is called to major allele only
transform_snp <- function(pat, binomial_n, levels = c(), get=c("levels", "transformed")){
    std_bases <- c("A", "C", "G", "T")
    dual_bases <- list("AC"="M",
                    "AG"="R",
                    "AT"="W",
                    "CG"="S",
                    "CT"="Y",
                    "GT"="K")
    pat <- toupper(pat)
    unique_bases <- unique(pat)
    pat_ambig_bases <- unique_bases[unique_bases %in% dual_bases]
    pat_std_bases <- sort(unique_bases[unique_bases %in% std_bases])

    ### Determine the transformation levels
    if (length(levels) == 0){
        if (length(unique_bases) > (binomial_n + 1))
            stop("Unable to transform: unique_bases > binomial_n")
        if (length(pat_std_bases) <= 1 && binomial_n == 1)
            warning("Unique_bases <= 1, N added as the second class")
            pat_std_bases <- c(pat_std_bases, "N")
        if (length(pat_ambig_bases) > 1)
            stop("Unable to transform: multiple ambiguity codes found")
        if (length(pat_std_bases) == 2 && length(pat_ambig_bases) == 1)
            stopifnot(dual_bases[paste0(pat_std_bases, collapse = "")] == pat_ambig_bases)
        if (length(pat_std_bases) == 1 && length(pat_ambig_bases) == 1){
            inferred_std_bases <- strsplit(names(dual_bases[which(dual_bases == pat_ambig_bases)]), split = "")[[1]]
            stopifnot(pat_std_bases %in% inferred_std_bases)
            pat_std_bases <- inferred_std_bases
        }
        if (binomial_n == 2 && length(pat_ambig_bases) == 0){
            pat_ambig_bases <- dual_bases[paste0(pat_std_bases, collapse = "")][[1]]
        }
        all_bases <- c(pat_std_bases[1], pat_ambig_bases, pat_std_bases[2])
        attr(all_bases, "order") <- seq_along(all_bases) - 1
        levels <- all_bases
    }
    if (all(get == "levels"))
        return(levels)
    ### Transformation steps
    temp_levels <- levels
    if (binomial_n == 1 && length(levels) > 2)
        temp_levels <- temp_levels[c(1,3)]
    stopifnot(length(temp_levels) == (binomial_n + 1))
    if (binomial_n == 2){
        pat <- gsub("N", levels[2], pat)
    }
    transformed_pat <- match(pat, temp_levels) - 1
    if (all(get == "transformed"))
        return(transformed_pat)
    return(list(transformed_pat = transformed_pat,
                levels = levels))
}
assignInNamespace("transform_snp", transform_snp, ns="minSNPs")

all_seqs <- read_fasta("simgsk.fasta")

metadata <- fread("metadata.csv")
colnames(metadata)[which(colnames(metadata) == "sample")] <- "isolate"
training_seqs <- all_seqs[metadata[partition == "training"]$isolate]
validation_seqs <- all_seqs[metadata[partition == "validation"]$isolate]

f_result_mcc_multi_country <- c("mcc_multi_country.tsv", "mcc_multi_country2.tsv")

f_result_simpson_group_country <- c("simpson_by_group_country.tsv", "simpson_by_group_country2.tsv")


result_m_c <- unlist(lapply(f_result_mcc_multi_country, process_result_file), recursive = FALSE)
result_s_c <- unlist(lapply(f_result_simpson_group_country, process_result_file), recursive = FALSE)


# from a list of results, incrementally combined together, 1, 1+2, 1+2+3, 1+2+3+4
result_s_c_p <- lapply(2:(length(result_s_c)), function(i){return(unique(unlist(result_s_c[1:i])))})
result_m_c_p <- lapply(2:(length(result_m_c)), function(i){return(unique(unlist(result_m_c[1:i])))})
set.seed(2024)
all_sampled_position <- sample(1:length(all_seqs[[1]]), 5*30)
split_sampled_position <- setNames(split(all_sampled_position, ceiling(seq_along(all_sampled_position)/5)), NULL)
random_sample <- lapply(2:length(split_sampled_position), function(i){return (unique(unlist(split_sampled_position[1:i])))})

to_iterate <- list("mcc_multi_country_it" = list(result = c(result_m_c, result_m_c_p),
                                target = "country"),
    "simpson_group_country_it" = list(result = c(result_s_c, result_s_c_p),
                                   target = "country"),
    "mcc_multi_country2region_it" = list(result = c(result_m_c, result_m_c_p),
                                      target = "region"),
    "simpson_group_country2region_it" = list(result = c(result_s_c, result_s_c_p),
                                          target = "region"),
    "random_country_it" = list(result = c(split_sampled_position, random_sample),
                               target = "country"),
    "random_region_it" = list(result = c(split_sampled_position, random_sample),
                                target = "region")
)

all_predicted_result <- list()
for (item_name in names(to_iterate)){
    item <- to_iterate[[item_name]]
    colnames(metadata)[which(colnames(metadata) == item$target)] <- "target"
    training_metadata <- metadata[partition == "training"]
    validation_metadata <- metadata[partition == "validation"]

    all_sets <- item$result
    training_results <- list()
    validation_results <- list()
    all_results <- list()
    raw_training_results <- list()
    raw_validation_results <- list()
    predicted_results <- list()
    i <- 1
    for (sset in all_sets){
        print(paste0("CURRENTLY RUNNING: ", item_name, " ", i))
        i <- i + 1
        # train balk based on training data
        temp_balk <- train_balk(training_seqs, sset, training_metadata)
        training_result <- as.character(predict_balk(temp_balk, type = "class"))

        # validate balk
        val_seq <- generate_pattern(validation_seqs, sset)
        validation_result <- as.character(predict_balk(temp_balk, val_seq, snp_id = sset, type = "class"))
        
        # calculate mcc
        training_r_with_truth <- data.table(isolate=names(training_seqs), predicted_target = training_result)
        training_r_with_truth <- merge(training_r_with_truth, training_metadata[,.(isolate, target)], by = "isolate")
        training_mcc <- mcc_calculation(training_r_with_truth)
        training_mcc_raw <- mcc_calculation(training_r_with_truth, is_multi = FALSE)

        validation_r_with_truth <- data.table(isolate=names(validation_seqs), predicted_target = validation_result)
        validation_r_with_truth <- merge(validation_r_with_truth, validation_metadata[,.(isolate, target)], by = "isolate")
        validation_mcc <- mcc_calculation(validation_r_with_truth)
        validation_mcc_raw <- mcc_calculation(validation_r_with_truth, is_multi = FALSE)

        all_r_with_truth <- rbind(training_r_with_truth, validation_r_with_truth)
        all_mcc <- mcc_calculation(all_r_with_truth)

        training_mcc$set <- paste(sset, collapse = "_")
        validation_mcc$set <- paste(sset, collapse = "_")
        all_mcc$set <- paste(sset, collapse = "_")
        training_mcc_raw$set <- paste(sset, collapse = "_")
        validation_mcc_raw$set <- paste(sset, collapse = "_")

        training_results[[paste(sset, collapse = "_")]] <- training_mcc
        validation_results[[paste(sset, collapse = "_")]] <- validation_mcc
        all_results[[paste(sset, collapse = "_")]] <- all_mcc
        raw_training_results[[paste(sset, collapse = "_")]] <- training_mcc_raw
        raw_validation_results[[paste(sset, collapse = "_")]] <- validation_mcc_raw

        all_r_with_truth$set <- paste(sset, collapse = "_")
        predicted_results[[paste(sset, collapse = "_")]] <- all_r_with_truth
    }
    temp <- rbindlist(predicted_results)
    temp$snp_source <- item_name
    all_predicted_result[[item_name]] <- temp
    training_results <- rbindlist(training_results)
    training_results$n_snps <- lengths(strsplit(training_results$set, split = "_"))
    validation_results <- rbindlist(validation_results)
    validation_results$n_snps <- lengths(strsplit(training_results$set, split = "_"))
    all_results <- rbindlist(all_results)
    all_results$n_snps <- lengths(strsplit(training_results$set, split = "_"))
    training_results$type <- "training"
    validation_results$type <- "validation"
    all_results$type <- "all"
    combined_results <- rbind(training_results, validation_results, all_results)
    raw_training_results <- rbindlist(raw_training_results)
    raw_validation_results <- rbindlist(raw_validation_results)

    fwrite(combined_results, paste0(item_name, "_combined_results_balk.csv"), row.names = FALSE)
    fwrite(raw_training_results, paste0(item_name, "_raw_training_results_balk.csv"), row.names = FALSE)
    fwrite(raw_validation_results, paste0(item_name, "_raw_validation_results_balk.csv"), row.names = FALSE)

    colnames(metadata)[which(colnames(metadata) == "target")] <- item$target
}


all_predicted_v_truth <- rbindlist(all_predicted_result)
all_predicted_v_truth$partition <- metadata$partition[match(all_predicted_v_truth$isolate, metadata$isolate)]
fwrite <- fwrite(all_predicted_v_truth, "all_predicted_v_truth_it.csv", row.names = FALSE)

all_combined_results_files <- list.files(pattern = "it_combined_results_balk.csv")
all_combined_results <- lapply(all_combined_results_files, function(f){
    res <- fread(f)
    res$snp_source <- gsub("_it_combined_results_balk.csv", "", f)
    return(res)
    })
all_combined_results <- rbindlist(all_combined_results)
all_combined_results[snp_source == "random_region"]$snp_source <- "random_country2region"
summary_result <- lapply(unique(all_combined_results$snp_source), function(source){
    if (grepl("region", source)){
        type <- "region"
    } else {
        type <- "country"
    }
    highest_set_1 <- all_combined_results[all_combined_results$snp_source == source & type == "training"][order(-mcc, n_snps)]$set[1]
    if (type == "country"){
        mcc_country <- all_combined_results[all_combined_results$snp_source == source & set == highest_set_1][match(type, c("training", "validation", "all"))]$mcc
        mcc_region <- all_combined_results[all_combined_results$snp_source == gsub("country", "country2region", source) & set == highest_set_1][match(type, c("training", "validation", "all"))]$mcc
    } else {
        mcc_region <- all_combined_results[all_combined_results$snp_source == source & set == highest_set_1][match(type, c("training", "validation", "all"))]$mcc
        mcc_country <- all_combined_results[all_combined_results$snp_source == gsub("country2region", "country", source) & set == highest_set_1][match(type, c("training", "validation", "all"))]$mcc
    }
    mcc_highest_set <- all_combined_results[all_combined_results$snp_source == source & set == highest_set_1][match(type, c("training", "validation", "all"))]$mcc
    data.frame(snp_source = source, highest_set = highest_set_1,
        training_mcc_country = mcc_country[1], validation_mcc_country = mcc_country[2], all_mcc_country = mcc_country[3],
        training_mcc_region = mcc_region[1], validation_mcc_region = mcc_region[2], all_mcc_region = mcc_region[3])
})
summary_result <- rbindlist(summary_result)
summary_result$n_snps <- lengths(strsplit(summary_result$highest_set, split = "_"))
fwrite(summary_result, "summarised_results_balk.csv", row.names = FALSE)


##### Generate heatmap
pairwise_mat <- function(data, is_region = FALSE){
    unique_columns <- unique(c(data$target, data$predicted_target))
    all_possible_columns <- c("Ethiopia", "Madagascar", "Sudan",
        "Afghanistan", "India", "Iran", "Bangladesh", "Bhutan",
        "Myanmar", "Thailand", "Cambodia", "Vietnam", "Malaysia",
        "Philippines", "Indonesia", "Papua New Guinea", "China",
        "Brazil", "Colombia", "Mexico", "Peru")
    if (is_region){
        all_possible_columns <- c("AF", "WAS", "WSEA", "ESEA", "MSEA", "OCE", "EAS", "LAM")
    }
    temp_col <- unique_columns[match(all_possible_columns, unique_columns)]
    columns <- temp_col[!is.na(temp_col)]
    result <- matrix(0, nrow = length(columns), ncol = length(columns))
    colnames(result) <- columns
    rownames(result) <- columns
    for (i in 1:length(columns)){
        for (j in 1:length(columns)){
            result[i, j] <- sum(data$target == columns[i] & data$predicted_target == columns[j])/sum(data$target == columns[i])
        }
    }
    return(result)
}

to_generate <- list("mcc_multi_single" = c(91113, 213179, 213507, 213300, 58712),
    "simpson_group_single" = c(221827, 183416, 117398, 118909, 213148),
    "mcc_multi_all" = tail(result_m_c_p, 1)[[1]],
    "simpson_group_all" = tail(result_s_c_p, 1)[[1]]
)

for (tg in names(to_generate)){
    setname <- paste0(to_generate[[tg]], collapse = "_")
    if(grepl("mcc", tg, fixed = TRUE)){
        snp_source <- c("mcc_multi_country_it", "mcc_multi_country2region_it")
    }else{
        snp_source <- c("simpson_group_country_it", "simpson_group_country2region_it")
    }
    for (ssource in snp_source){
        # ALL partition
        country_or_region <- ifelse(grepl("region", ssource), "region", "country")
        mcc_or_simpson <- ifelse(grepl("mcc", tg, fixed = TRUE), "mcc", "simpson")
        single_snp_or_all <- ifelse(grepl("single", tg, fixed = TRUE), "single", "all")
        print(paste0(mcc_or_simpson, "_", country_or_region, "_", single_snp_or_all, "_all.png")) #png(paste0(mcc_or_simpson, "_", country_or_region, "_", single_snp_or_all, "_all.png"))
        temp_dt <- all_predicted_v_truth[set == setname & snp_source == ssource][,.(target, predicted_target)]
        mat <- pairwise_mat(temp_dt, is_region = grepl("region", ssource))
        heat_col <- rev(heat.colors(100))
        heatmap(mat, Rowv = NA, Colv = NA, col=heat_col, margins = c(20, 20),
            main = "Heatmap of pairwise comparison of predicted vs truth", xlab = "", ylab = "")
        mtext("Truth", side = 4, line = -9)
        mtext("Predicted", side = 1, line = -2)
        legend(x="right", legend=c("0", "0.25", "0.5", "0.75", "1.0"),fill=heat_col[c(1,25,50,75,100)])
        #dev.off()
        invisible(readline(prompt="Press [enter] to continue"))
        # Validation partition
        print(paste0(mcc_or_simpson, "_", country_or_region, "_", single_snp_or_all, "_validation.png"))#png(paste0(mcc_or_simpson, "_", country_or_region, "_", single_snp_or_all, "_all.png"))
        temp_dt <- all_predicted_v_truth[set == setname & snp_source == ssource][partition == "validation"][,.(target, predicted_target)]
        mat <- pairwise_mat(temp_dt, is_region = grepl("region", ssource))
        heat_col <- rev(heat.colors(100))
        heatmap(mat, Rowv = NA, Colv = NA, col=heat_col, margins = c(20, 20),
            main = "Heatmap of pairwise comparison of predicted vs truth", xlab = "", ylab = "")
        mtext("Truth", side = 4, line = -9)
        mtext("Predicted", side = 1, line = -2)
        legend(x="right", legend=c("0", "0.25", "0.5", "0.75", "1.0"),fill=heat_col[c(1,25,50,75,100)])
        #dev.off()
        invisible(readline(prompt="Press [enter] to continue"))
    }
}



metadata <- fread("metadata.csv")
data <- fread("all_predicted_v_truth.csv")
data$partition <- metadata$partition[match(data$isolate, metadata$sample)]
longest_set <- unique(data[,.(snp_source, set = set)])[, .SD[which.max(nchar(set)) ], by = snp_source]
i <- 4
print(longest_set[i]$snp_source)
temp_dt <- data[snp_source == longest_set[i]$snp_source & set == longest_set[i]$set][,.(target, predicted_target)]#[partition == "validation"][,.(target, predicted_target)]
mat <- pairwise_mat(temp_dt,
    is_region = grepl("region", longest_set[i]$snp_source))
heat_col <- rev(heat.colors(100))
heatmap(mat, Rowv = NA, Colv = NA, col=heat_col, margins = c(10, 10),
    main = "Heatmap of pairwise comparison of predicted vs truth", xlab = "", ylab = "")
mtext("Truth", side = 4, line = -9)
mtext("Predicted", side = 1, line = 5)
legend(x="right", legend=c("0", "0.25", "0.5", "0.75", "1.0"),fill=heat_col[c(1,25,50,75,100)])

"mcc_multi_country_single"
c(91113, 213179, 213507, 213300, 58712)
"simpson_group_country_single" c(221827, 183416, 117398, 118909, 213148)

longest_set <- unique(data[,.(snp_source, set = set)])[set %in% c(
    paste0(c(91113, 213179, 213507, 213300, 58712), collapse = "_"),
    paste0(c(221827, 183416, 117398, 118909, 213148), collapse = "_")
)]