MAPK1 QSAR Modeling and Argemone mexicana Phytochemical Virtual Screening
Overview
This repository contains the Python workflow used for a quantitative structure–activity relationship (QSAR) study of MAPK1 (CHEMBL4040), followed by virtual screening of phytochemicals reported for Argemone mexicana.
The workflow retrieves activity records from the ChEMBL web service, curates IC50 measurements, calculates Morgan fingerprints, trains a Random Forest regression model to predict pIC50, evaluates model performance, generates an applicability-domain plot, and predicts pIC50 values for the phytochemical input file.
> **Important:** This README describes the workflow represented in the submitted script. Before using this repository as the final publication version, resolve the code issues listed under “Checks required before release” and rerun the workflow to confirm that the published code reproduces the reported results.
Workflow
Retrieve activity records: Query the ChEMBL activity API for target `CHEMBL4040` and retrieve paginated activity records.
Curate activity data: Retain records with `standard_type = IC50`, `standard_units = nM`, and `standard_relation = =`. Convert the activity values to numeric values, discard non-positive IC50 values and records without canonical SMILES, and calculate pIC50 as:
[
\mathrm{pIC}{50}=9-\log{10}(\mathrm{IC}_{50}\text{ in nM})
]
Aggregate repeated measurements: Group records by ChEMBL molecule ID and calculate median IC50, median pIC50, and the number of measurements per compound.
Generate molecular descriptors: Calculate radius-2 Morgan bit fingerprints with 2,048 bits using RDKit.
Split the dataset: Use an 80:20 training/test split with `random_state=42`.
Train the QSAR model: Fit a `RandomForestRegressor` with 500 trees, `max_features="sqrt"`, and `random_state=42`.
Evaluate predictions: Calculate test-set R², RMSE, and MAE; perform shuffled five-fold cross-validation; and generate an actual-versus-predicted plot.
Assess applicability domain: Calculate maximum Tanimoto similarity of test compounds to training compounds and plot similarity against absolute prediction error.
Screen phytochemicals: Read the phytochemical spreadsheet, generate the same type of Morgan fingerprints, predict pIC50 values, and rank the compounds by predicted pIC50.
Repository contents
Suggested organization:
```text
.
├── README.md
├── qsar_mapk1.py
├── requirements.txt
└── .gitignore
```

Requirements
The script uses Python and the following packages:
`numpy`
`pandas`
`matplotlib`
`seaborn`
`requests`
`rdkit`
`scikit-learn`
`openpyxl` (for Excel input/output)
The workflow also requires an internet connection to query the ChEMBL API. Record and publish the exact Python and package versions used for the final analysis after verifying the environment in which the workflow runs.
Example installation with Conda (recommended for RDKit):
```bash
conda create -n mapk1-qsar python=3.11
conda activate mapk1-qsar
conda install -c conda-forge numpy pandas matplotlib seaborn requests rdkit scikit-learn openpyxl
```
Test the environment and adjust package versions as necessary before reporting it as the validated environment for this study.
Input data
ChEMBL activity data
The script queries the ChEMBL activity API for target ID `CHEMBL4040`. The activity records are retrieved at runtime, so the resulting dataset may change if ChEMBL records are updated. For reproducibility, retain the curated dataset used for the final manuscript, record the retrieval date, and document any additional inclusion/exclusion criteria.
Argemone mexicana phytochemical data
The script expects an Excel file named `Argemone mexicana.xlsx` with a column containing compound names (`Phytochemical name`) and a SMILES column (`SMILES` or `Smile`). Place the file in the configured input-data directory and confirm that the column names match the script.

Running the workflow
Install the required packages.
Place the QSAR script and phytochemical spreadsheet in the intended project folders.
Update the hard-coded Windows output path in the script to a project-relative path or another path suitable for your computer.
Confirm that the input spreadsheet contains the expected columns.
Run the script:
```bash
   python qsar_mapk1.py
   ```
Check the terminal output and generated tables/figures. Confirm that all outputs are produced and that the reported values match the final manuscript.
The script currently writes Excel files and plots to the output directory configured inside the code. It should be updated to use the repository's `data/` and `results/` directories before release.
Model evaluation
The script reports the following metrics:
R²: coefficient of determination on the held-out test set and across cross-validation folds.
RMSE: root mean squared error.
MAE: mean absolute error.
Applicability-domain analysis: maximum Tanimoto similarity to training compounds plotted against absolute prediction error.
The script currently labels the mean of the five fold-wise R² values as “Q²”. State the exact definition used in the manuscript and verify that it matches the calculation. A mean of fold-wise R² values is not necessarily identical to a pooled cross-validated Q² calculation.
