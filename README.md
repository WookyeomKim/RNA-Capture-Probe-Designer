# Hybridization Capture Probe Designer for RNA Extraction

This repository provides a Python tool for designing sequence-specific hybridization capture probes. It is particularly effective for designing a minimal number of DNA capture probes, such as one or two, for direct RNA extraction from plasma and whole blood, as described in our research.

The algorithm performs a single-pass screening of genomic sequences to rank candidate probes based on mismatch tolerance, GC content, and thermodynamic binding stability calculated using NUPACK. Accessibility, defined as the probability that contiguous 10 nt within the capture site are fully unpaired, is also calculated using the RNAplfold module in ViennaRNA.

A detailed description of the algorithm, its underlying rationale, and experimental validation will be provided in our manuscript, currently under submission.


## 1. Key Features
- **Single-Pass Ranking:** Evaluates all potential sites once and ranks them by performance.
- **Memory Optimization:** Efficiently handles large FASTA datasets using batch processing and incremental history saving.
- **Thermodynamic Validation:** Integrates NUPACK for precise binding energy (Delta G) calculations.
- **Accessibility Scan:** Uses `RNAplfold` to ensure probes target accessible regions of the target RNA.

## 2. Prerequisites
Before running the script, ensure the following are installed on your system:

### **External Software**
- **ViennaRNA Package (RNAplfold):** This is required for accessibility scanning. 
  - Install via Conda: `conda install -c bioconda viennarna`
  - Ensure `RNAplfold` is in your system's PATH.

### **Python Packages**
Install the necessary Python libraries using the provided `requirements.txt`:
```bash
pip install -r requirements.txt