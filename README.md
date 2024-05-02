
# Drug Repositioning Analysis Tool

This repository contains a Python script for analyzing drug repositioning opportunities using gene expression data. The script processes specific data formats and integrates several data sources to identify potential new uses for existing drugs based on their impact on gene expression in various cell lines.

## Features

- Filtering and processing of signature and compound information.
- Batch processing of gene expression data through `.gctx` files.
- Calculation of XSum scores for differential gene expression.
- Integration of compound aliases and clinical phase data for potential repurposing.
- Output of repositioned drug candidates with relevant metadata.

## Installation

To run the script, you need to install the required Python packages. You can install them using the following command:

```bash
pip install -r requirements.txt
```

## Usage

The script uses command-line arguments for configuration. Here’s how to run the script:

```bash
python drug_repositioning.py --sig_info PATH_TO_SIGINFO --gene_info PATH_TO_GENEINFO --query PATH_TO_QUERY --batch_size BATCH_SIZE --trt_cp PATH_TO_GCTX --Xsum_num TOP_N --cmp_info PATH_TO_COMPOUND_INFO --drug_info PATH_TO_DRUG_INFO --output OUTPUT_FILENAME --cell CELL_LINE
```

### Arguments

- `--sig_info`: Path to the signature information file.
- `--gene_info`: Path to the gene information file.
- `--query`: Path to the Excel file containing the query genes.
- `--batch_size`: Number of signatures to process in each batch.
- `--trt_cp`: Path to the `.gctx` file for treatment compound data.
- `--Xsum_num`: Number of top and bottom expressions to consider for XSum scoring.
- `--cmp_info`: Path to the compound information file.
- `--drug_info`: Path to the drug information file containing repurposing drugs.
- `--output`: Filename for the output CSV containing repositioned drugs.
- `--cell`: Cell line to be analyzed.

## Contributing

Contributions to this project are welcome. You can contribute by:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Becoming a maintainer

## License

Distributed under the MIT License. See `LICENSE` for more information.

## Contact

Your Name – Your Email

Project Link: https://github.com/yourusername/yourrepositoryname

