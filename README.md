# minSNPs_pv4

## Deriving SNP sets to discriminate _P. vivax_ samples from east southeast asia, and Vietnam with "mcc" and "percent"
- Scripts used:
    - [ESEA samples](scripts/)
    - [Vietnam samples](scripts/)
- Results:


## Deriving SNP sets with generalised discriminatory power with "simpson" for comparison
- Scripts used:

- Results:

## Deriving SNP sets to discriminate _P. vivax_ according to region of infection with "mcc_multi" and "simpson_by_group"
- Scripts used:

- Results:

## Deriving SNP sets to discriminate _P. vivax_ according to country of infection with "mcc_multi" and "simpson_by_group"
- Scripts used:

- Results:

## Combined analysis:
- Scripts used:

- Results:

## Using the "mcc_multi" and "simpson_by_group" SNPs to build BALK classifier
- Scripts used:
    - `scripts/BALK/`

- Results:
    - Heatmaps are in `results/BALK/heatmaps/`, in the files with following format:
        - <mcc_multi | simpson_group>\_<country | country2region>\_<all | val>(_single).png
    - The final classifiers is in `results/BALK/`