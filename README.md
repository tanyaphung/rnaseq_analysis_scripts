# rnaseq_analysis_scripts
This repo contains scripts relating to analysis of RNAseq data

## Measuring expression level from RNAseq data
- Scripts written in Python following: http://training.scicomp.jic.ac.uk/docs/hpc_rnaseq_course_book/expression.html
- `measure_expression_level.py`
```
python measure_expression_level.py -h
usage: measure_expression_level.py [-h] --featureCounts_file
                                   FEATURECOUNTS_FILE --outfile OUTFILE
                                   --normalization NORMALIZATION

Compute TPM (transcripts per million) from featureCounts output

optional arguments:
  -h, --help            show this help message and exit
  --featureCounts_file FEATURECOUNTS_FILE
                        Input the path to the featureCounts output
  --outfile OUTFILE     Input the path to the output file after computing TPM
  --normalization NORMALIZATION
                        Input the types of normalization. Currently only
                        support either tpm or fpkm
```
