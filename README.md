# minSNPs_pv4

*note:*
The SNPs from the different chromosomes are concatenated into a single file. The selected SNPs are the index in the concatenated FASTA file. `ref_pv4.tsv` can be used to find the chromosome and the relative position of the SNP in the chromosome. 

## Deriving SNP sets to discriminate _P. vivax_ samples from east southeast asia, and Vietnam with "mcc" and "percent"
- Results:
    - Vietnam: `results/minSNPs_output/percent_vietnam.tsv`, `results/minSNPs_output/mcc_vietnam.tsv`
    - ESEA: `results/minSNPs_output/percent_esea.tsv`, `results/minSNPs_output/mcc_esea.tsv`

## Deriving SNP sets with generalised discriminatory power with "simpson" for comparison
- Results: `results/minSNPs_output/simpson_old.tsv`

## Deriving SNP sets to discriminate _P. vivax_ according to region and country of infection with "mcc_multi" and "simpson_by_group"
- Results: 
    - Region of infection:
        - `results/minSNPs_output/mcc_multi_region.tsv`
        - `results/minSNPs_output/simpson_by_group_region.tsv` and `results/minSNPs_output/simpson_by_group_region2.tsv`
    - Country of infection:
        - `results/minSNPs_output/mcc_multi_country.tsv` and `results/minSNPs_output/mcc_multi_country2.tsv`
        - `results/minSNPs_output/simpson_by_group_country.tsv` and `results/minSNPs_output/simpson_by_group_country2.tsv`

## Combined analysis:
- Scripts used: `scripts/minSNPs/inf_summ.r`
- Results:
    - `results\combined_analysis\mcc_v_percent_esea_vn.csv`
    - `results\combined_analysis\mcc_v_percent_esea_vn_raw.csv`
    - `results\combined_analysis\mcc_simpson_by_group_v_simpson_result.csv`
    - `results\combined_analysis\mcc_simpson_by_group_v_simpson_result_raw_region.csv`
    - `results\combined_analysis\mcc_simpson_by_group_v_simpson_result_raw_country.csv`

## Using the "mcc_multi" and "simpson_by_group" SNPs to build BALK classifier
- Scripts used:
    - `scripts/BALK/BALK_experiment.r`

- Results:
    - The results to compare model built with SNPs selected mcc_multi, simpson_by_group and random are in `results/BALK/`
    - `results\BALK\all_predicted_v_truth_it.csv`
    - `results\BALK\mcc_multi_country_it_combined_results_balk.csv`
    - `results\BALK\mcc_multi_country_it_raw_training_results_balk.csv`
    - `results\BALK\mcc_multi_country_it_raw_validation_results_balk.csv`
    - `results\BALK\mcc_multi_country2region_it_combined_results_balk.csv`
    - `results\BALK\mcc_multi_country2region_it_raw_training_results_balk.csv`
    - `results\BALK\mcc_multi_country2region_it_raw_validation_results_balk.csv`
    - `results\BALK\random_country_it_combined_results_balk.csv`
    - `results\BALK\random_country_it_raw_training_results_balk.csv`
    - `results\BALK\random_country_it_raw_validation_results_balk.csv`
    - `results\BALK\random_region_it_combined_results_balk.csv`
    - `results\BALK\random_region_it_raw_training_results_balk.csv`
    - `results\BALK\random_region_it_raw_validation_results_balk.csv`
    - `results\BALK\simpson_group_country_it_combined_results_balk.csv`
    - `results\BALK\simpson_group_country_it_raw_training_results_balk.csv`
    - `results\BALK\simpson_group_country_it_raw_validation_results_balk.csv`
    - `results\BALK\simpson_group_country2region_it_combined_results_balk.csv`
    - `results\BALK\simpson_group_country2region_it_raw_training_results_balk.csv`
    - `results\BALK\simpson_group_country2region_it_raw_validation_results_balk.csv`
    - `results\BALK\summarised_results_balk.csv`
    - Heatmaps are in `results/BALK/heatmaps/`, in the files with following format:
        - <mcc_multi | simpson_group>\_<country | country2region>\_<all | val>(_single).png