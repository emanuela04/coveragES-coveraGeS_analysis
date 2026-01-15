# Low-Coverage Regions Analysis in NGS data

This directory contains the scripts developed and used for the analysis of **low-coverage regions (LCRs)** in **paired whole-exome sequencing (ES) and whole-genome sequencing (GS)** data generated from the same individuals.

The primary goal of this study is to assess the **distribution, recurrence, and concordance** of low-coverage regions between ES and GS obtained using [Mosdepth](https://github.com/brentp/mosdepth), and to evaluate their **potential impact on clinically relevant genes** with respect to diagnostic interpretation.

```
mosdepth \
  --by gencode.v35.protein_coding_CDS.bed \
  --thresholds 1,10,20,30 \
  --fasta GRCh38_full_analysis_set_plus_decoy_hla.fa \
  sample \
  sample.cram
```


Generate coding intervals (CDS) from Gencode v35 GTF: Protein-coding CDS intervals were extracted from the Gencode v35 GTF and converted to BED (0-based start, 1-based end).


```
zgrep -v "^##" gencode.v35.annotation.gtf.gz \
  | awk -F'\t' '
      $3=="CDS" && $0 ~ /transcript_type "protein_coding"/ {
        # BED format: chr, start(0-based), end, attributes
        print $1, $4-1, $5, $9
      }
    ' OFS='\t' \
  > gencode.v35.protein_coding_CDS.bed
```


A comparative analysis between ES and GS was performed on **14 paired samples**, allowing a direct evaluation of technology-specific versus intrinsic low-coverage regions.

In addition, a dedicated analysis was conducted to investigate **batch effects** in low-coverage profiles using **four independent ES sequencing batches**, in order to characterize systematic coverage biases and batch-specific behaviors.

Clinical annotation of low-coverage genomic intervals was further enriched using an updated **uncoverappLib-based annotation workflow**, enabling the identification of **potentially clinical relevant variants** located within low-coverage regions.  
The updated version of [uncoverappLib](https://github.com/emanuela04/uncoverappLib/tree/v2.1.0) overcomes the performance limitations of the previous implementation while preserving all its core functionalities.

Overall, this repository provides a reproducible framework for the systematic characterization of low-coverage regions and their relevance in clinical sequencing diagnostics.

