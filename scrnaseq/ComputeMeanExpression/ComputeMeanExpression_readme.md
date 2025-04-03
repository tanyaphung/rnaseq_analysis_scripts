## Purpose: 
- To compute mean (or specificity) per gene from scRNAseq data
- Demo with the data from Velmeshev et al. 2023 (https://cellxgene.cziscience.com/collections/bacccb91-066d-4453-b70e-59de0b4598cd)

## Step 1: Create a conda environment
- yml file is provided in `resources/scrna_env.yml`
```
conda env create --name scrna --file=resources/scrna_env.yml
```

## Step 2: Download scRNAseq data in h5ad format
- Go to https://cellxgene.cziscience.com/collections/bacccb91-066d-4453-b70e-59de0b4598cd
- Click on `Download` and `.h5ad`
- Download with `wget`: 
```
wget https://datasets.cellxgene.cziscience.com/8110794f-d8c7-4611-8d38-67f5e22e41fb.h5ad
```

## Step 3: Explore the data
- If it is the first time I am analyzing the data, I would like to first explore the data to know:
    - what kinds of information are available in the `adata.obs`? 
    - is that count raw or processed? 
    - how many genes? 
- To do this, run the script `explore.py`
```
python explore.py -f 8110794f-d8c7-4611-8d38-67f5e22e41fb.h5ad -o explore_out.md
```

- Examining the output `explore_out.md`: 
    - we see that there are several brain regions (column `tissue`):
    ```
    frontal cortex
    cingulate cortex
    insular cortex
    cerebral cortex
    temporal cortex
    lateral ganglionic eminence
    caudal ganglionic eminence
    medial ganglionic eminence
    ganglionic eminence
    primary motor cortex
    ```
    - there are several developmental stages (column `development_stage`)
    ```
    12-year-old stage
    14-year-old stage
    22-year-old stage
    21-year-old stage
    19-year-old stage
    13-year-old stage
    15-year-old stage
    8-year-old stage
    6-year-old stage
    4-year-old stage
    39-year-old stage
    54-year-old stage
    34-year-old stage
    44-year-old stage
    2-year-old stage
    infant stage
    ninth LMP month stage
    2-month-old stage
    3-year-old stage
    1-month-old stage
    1-year-old stage
    fifth LMP month stage
    eighth LMP month stage
    sixth LMP month stage
    seventh LMP month stage
    5-year-old stage
    3-month-old stage
    16-year-old stage
    fourth LMP month stage
    45-year-old stage
    53-year-old stage
    28-year-old stage
    newborn stage (0-28 days)
    4-month-old stage
    5-month-old stage
    10-year-old stage
    20-year-old stage
    40-year-old stage
    17-year-old stage
    25-year-old stage
    ```

## Step 4: Subset the data to each brain region and each developmental stage
- It could be the case that not all developmental stage is available for all brain region. 
- I would loop through each tissue and each developmental age to subset. If the file is empty, then remove them
- Example: 
```
for region in FrontalCortex CingulateCortex InsularCortex CerebralCortex TemporalCortex LateralGanglionicEminence CaudalGanglionicEminence MedialGanglionicEminence GanglionicEminence PrimaryMotorCortex
do
for age in 12Y 14Y 22Y 21Y 19Y 13Y 15Y 8Y 6Y 4Y 39Y 54Y 34Y 44Y 2Y Infant 9LMPMonth 2Month 3Y 1Month 1Y 5LMPMonth 8LMPMonth 6LMPMonth 7LMPMonth 5Y 3Month 16Y 4LMPMonth 45Y 53Y 28Y Newborn 4Month 5Month 10Y 20Y 40Y 17Y 25Y
do
mkdir ${region}_${age}
python subset.py -d {directory} -f 8110794f-d8c7-4611-8d38-67f5e22e41fb.h5ad -id ${region}_${age}
done
done
```

- After this step, one would have subset the data for per region and per developmental age (if there are no cells in the subsetted data, remove the ID from further analyses). 
- **NOTES**: the script `subset.py` requires `single_cell_helper_functions_v3.py` so make sure that both scripts are in the same directory. 

## Step 5: compute mean per gene
- To compute mean and specificity per gene, TODO: edit the script `https://github.com/tanyaphung/scrnaseq_viewer/blob/main/code/postprocessing/compute_sumstat_magma.py` to remove the `per cell type`. How to run this script: 
    ```
    python compute_sumstat_magma.py --h5ad {input} --outdir {params.outdir} --ct_colname {params.ct_colname}
    ```