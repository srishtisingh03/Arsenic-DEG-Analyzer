# Arsenic-DEG-Analyzer

A Python-based bioinformatics tool for analyzing differential gene expression (DEG) data from arsenic-stress transcriptomics studies.

## Features

- Detects upregulated genes
- Detects downregulated genes
- Generates DEG summary statistics
- Creates volcano plots
- Exports DEG reports

## Requirements

- Python 3.x
- pandas
- matplotlib
- numpy

Install dependencies:

pip install -r requirements.txt

## Input Format

Example CSV:

Gene,logFC,pvalue
OsNIP2;1,2.5,0.001
OsABCC1,-1.8,0.02
OsPCS1,3.1,0.0001

## Usage

Run:

python deg_analyzer.py

## Output Files

The tool generates:

- upregulated_genes.csv
- downregulated_genes.csv
- summary.txt
- volcano_plot.png

## Applications

- Rice arsenic transcriptomics
- RNA-seq DEG analysis
- Candidate gene identification
- Bioinformatics workflows

## Author

Srishti Singh
