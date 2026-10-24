# ==========================================
# Cell 1: Environment Setup & Imports
# ==========================================
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.base import clone
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, train_test_split
# Define directory path and ensure it exists
output_dir = r"C:\Users\Admin\Argemone Mexicana"
os.makedirs(output_dir, exist_ok=True)
print("Python:", sys.version)
print("Environment: ARMEONA-QSAR")
print(f"Target Directory: {output_dir}")
print("Ready for ChEMBL-ARMEONA QSAR & Validation")
# ==========================================
# Cell 2: API Initial Connection Test
# ==========================================
target_id = "CHEMBL4040"  # MAPK1 target ID
url = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
params = {"target_chembl_id": target_id, "limit": 1000}
response = requests.get(url, params=params, timeout=60)
response.raise_for_status()
data = response.json()
activities = pd.DataFrame(data["activities"])
print("Records retrieved:", len(activities))
print("Columns:", len(activities.columns))
print("Target:", target_id)
# ==========================================
# Cell 3: Fetching All Target Activity Records
# ==========================================
all_activities = []
offset = 0
limit = 1000
while True:
    params = {"target_chembl_id": target_id, "limit": limit, "offset": offset}
    r = requests.get(
        "https://www.ebi.ac.uk/chembl/api/data/activity.json",
        params=params,
        timeout=60,
    )
    r.raise_for_status()
    page = r.json()["activities"]
    all_activities.extend(page)
    print(f"Retrieved: {len(all_activities)}")
    if len(page) < limit:
        break
    offset += limit
target_raw = pd.DataFrame(all_activities)
print(f"\nTotal activity records for {target_id}:", len(target_raw))
raw_output_path = os.path.join(output_dir, f"{target_id}_Raw_Activities_Summary.xlsx")
target_raw.head(1000).to_excel(raw_output_path, index=False)
# ==========================================
# Cell 4: Inspecting Activity Types & Units
# ==========================================
print("Activity types:")
print(target_raw["standard_type"].value_counts(dropna=False))
print("\nStandard units:")
print(target_raw["standard_units"].value_counts(dropna=False))
# ==========================================
# Cell 5: Filtering and Cleaning IC50 Data
# ==========================================
target_ic50 = target_raw[
    (target_raw["standard_type"] == "IC50")
    & (target_raw["standard_units"] == "nM")
    & (target_raw["standard_relation"] == "=")
].copy()
target_ic50["IC50_nM"] = pd.to_numeric(
    target_ic50["standard_value"], errors="coerce"
)
target_ic50 = target_ic50[
    (target_ic50["IC50_nM"] > 0) & (target_ic50["canonical_smiles"].notna())
].copy()
target_ic50["pIC50"] = 9 - np.log10(target_ic50["IC50_nM"])
print("Clean IC50 records:", len(target_ic50))
print("Unique compounds:", target_ic50["molecule_chembl_id"].nunique())
# ==========================================
# Cell 6: Grouping Duplicates & Calculating Medians
# ==========================================
target_final = (
    target_ic50.groupby("molecule_chembl_id")
    .agg(
        canonical_smiles=("canonical_smiles", "first"),
        median_IC50_nM=("IC50_nM", "median"),
        median_pIC50=("pIC50", "median"),
        n_measurements=("pIC50", "count"),
    )
    .reset_index()
)
print("Final unique compounds:", len(target_final))
unique_output_path = os.path.join(output_dir, f"{target_id}_Unique_Compounds.xlsx")
target_final.to_excel(unique_output_path, index=False)
# ==========================================
# Cell 7: Generating Morgan Fingerprints
# ==========================================
fps = []
valid_rows = []
for i, smi in enumerate(target_final["canonical_smiles"]):
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        fps.append(np.array(fp))
        valid_rows.append(i)
X = np.array(fps)
y = target_final.loc[valid_rows, "median_pIC50"].to_numpy()
print("Valid molecules:", len(valid_rows))
print("Fingerprint shape:", X.shape)
print("Target shape:", y.shape)
# ==========================================
# Cell 8: Train/Test Splitting
# ==========================================
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.20, random_state=42
)
print("Training compounds:", len(X_train))
print("Test compounds:", len(X_test))
# ==========================================
# Cell 9: Training the Random Forest Model
# ==========================================
rf_model = RandomForestRegressor(
    n_estimators=500, max_features="sqrt", random_state=42, n_jobs=-1
)
rf_model.fit(X_train, y_train)
print("Random Forest model trained successfully.")
# ==========================================
# Cell 10: External Test-Set Evaluation & Plot
# ==========================================
y_pred_train = rf_model.predict(X_train)
y_pred_test = rf_model.predict(X_test)
test_r2 = r2_score(y_test, y_pred_test)
test_rmse = np.sqrt(mean_squared_error(y_test, y_pred_test))
test_mae = mean_absolute_error(y_test, y_pred_test)
print("\n--- External Test-Set Performance ---")
print(f"Test R²  : {test_r2:.3f}")
print(f"Test RMSE: {test_rmse:.3f}")
print(f"Test MAE : {test_mae:.3f}")
test_metrics_df = pd.DataFrame(
    {"Metric": ["R2 Score", "RMSE", "MAE"], "Value": [test_r2, test_rmse, test_mae]}
)
test_metrics_df.to_excel(os.path.join(output_dir, f"{target_id}_Test_Set_Metrics.xlsx"), index=False)
sns.set_theme(style="whitegrid")
plt.figure(figsize=(9, 8), dpi=300)
plt.scatter(y_train, y_pred_train, color="blue", alpha=0.7, edgecolors="k", linewidths=0.3, label="Training")
plt.scatter(y_test, y_pred_test, color="red", alpha=0.85, edgecolors="k", linewidths=0.3, label="Validation")
all_actuals = np.concatenate([y_train, y_test])
all_preds = np.concatenate([y_pred_train, y_pred_test])
min_val = min(np.min(all_actuals), np.min(all_preds)) - 0.2
max_val = max(np.max(all_actuals), np.max(all_preds)) + 0.2
plt.plot([min_val, max_val], [min_val, max_val], color="black", linewidth=1.8, linestyle="-", label="Perfect Prediction")
plt.xlim(min_val, max_val)
plt.ylim(min_val, max_val)
plt.xlabel("Actual pIC50", fontsize=12, labelpad=10)
plt.ylabel("Predicted pIC50", fontsize=12, labelpad=10)
plt.title(f"{target_id} QSAR Model: Actual vs Predicted pIC50 Values", fontsize=14, pad=12, fontweight="medium")
plt.legend(loc="upper left", frameon=True, fontsize=10)
plt.tight_layout()
plot_output_path = os.path.join(output_dir, f"{target_id}_Actual_vs_Predicted_Train_Validation_Plot.png")
plt.savefig(plot_output_path, dpi=300)
plt.show()
# ==========================================
# Cell 11: 5-Fold Cross-Validation & Q² Calculation
# ==========================================
kf = KFold(n_splits=5, shuffle=True, random_state=42)
cv_r2_scores = []
cv_rmse_scores = []
cv_mae_scores = []
for train_idx, val_idx in kf.split(X):
    model = clone(rf_model)
    model.fit(X[train_idx], y[train_idx])
    pred = model.predict(X[val_idx])
    cv_r2_scores.append(r2_score(y[val_idx], pred))
    cv_rmse_scores.append(np.sqrt(mean_squared_error(y[val_idx], pred)))
    cv_mae_scores.append(mean_absolute_error(y[val_idx], pred))
q2_score = np.mean(cv_r2_scores)
cv_rmse_mean = np.mean(cv_rmse_scores)
cv_mae_mean = np.mean(cv_mae_scores)
print("\n--- 5-Fold Cross-Validation Metrics ---")
print(f"Cross-Validation Q² (Mean R²): {q2_score:.3f}")
print(f"Mean CV RMSE                : {cv_rmse_mean:.3f}")
print(f"Mean CV MAE                 : {cv_mae_mean:.3f}")
cv_results_df = pd.DataFrame({
    "Fold": [1, 2, 3, 4, 5],
    "R2_Score (CV)": cv_r2_scores,
    "RMSE": cv_rmse_scores,
    "MAE": cv_mae_scores
})
cv_results_df.to_excel(os.path.join(output_dir, f"{target_id}_5Fold_CV_Results.xlsx"), index=False)
# ==========================================
# Cell 12: Comprehensive Model Validation Summary
# ==========================================
validation_summary = pd.DataFrame({
    "Validation Parameter": ["External Test R²", "Cross-Validation Q²", "Test RMSE", "CV RMSE", "Test MAE", "CV MAE"],
    "Value": [test_r2, q2_score, test_rmse, cv_rmse_mean, test_mae, cv_mae_mean]
})
print("\n==========================================")
print("FINAL MODEL VALIDATION PARAMETERS SUMMARY:")
print("==========================================")
print(validation_summary.to_string(index=False))
validation_summary.to_excel(os.path.join(output_dir, f"{target_id}_Model_Validation_Summary.xlsx"), index=False)
# ==========================================
# Cell 13: Applicability Domain Scatter Plot Generation (No True/False Flags)
# ==========================================
train_indices = train_test_split(
    np.arange(len(target_final)), test_size=0.20, random_state=42
)[0]
test_indices = train_test_split(
    np.arange(len(target_final)), test_size=0.20, random_state=42
)[1]
train_count_fps = [
    rdMolDescriptors.GetMorganFingerprint(Chem.MolFromSmiles(smi), radius=2)
    for smi in target_final.iloc[train_indices]["canonical_smiles"]
    if Chem.MolFromSmiles(smi) is not None
]
test_count_fps = []
for idx in test_indices:
    mol = Chem.MolFromSmiles(target_final.iloc[idx]["canonical_smiles"])
    if mol is not None:
        test_count_fps.append(rdMolDescriptors.GetMorganFingerprint(mol, radius=2))
max_similarities = []
for test_fp in test_count_fps:
    similarities = DataStructs.BulkTanimotoSimilarity(test_fp, train_count_fps)
    max_similarities.append(max(similarities))
absolute_errors = np.abs(y_test - y_pred_test)
ad_results = pd.DataFrame({
    "Max_Tanimoto_Similarity": max_similarities,
    "Absolute_Error": absolute_errors,
})
ad_threshold = np.mean(max_similarities) - np.std(max_similarities)
print("\n--- Applicability Domain (AD) Plot Generation ---")
print(f"Mean Max Tanimoto Similarity: {np.mean(max_similarities):.3f}")
print(f"AD Threshold (Similarity Cutoff): {ad_threshold:.3f}")
# Plot AD Tanimoto Similarity vs Absolute Error cleanly without True/False hue labels
plt.figure(figsize=(7, 5), dpi=300)
sns.scatterplot(
    data=ad_results,
    x="Max_Tanimoto_Similarity",
    y="Absolute_Error",
    color="purple",
    alpha=0.6,
    edgecolor="k",
    linewidth=0.3
)
plt.axvline(
    x=ad_threshold,
    color="red",
    linestyle="--",
    linewidth=1.5,
    label=f"AD Threshold ({ad_threshold:.2f})",
)
plt.title(
    f"Applicability Domain: Tanimoto Similarity vs Absolute Error ({target_id})",
    fontsize=12,
)
plt.xlabel("Max Tanimoto Similarity to Training Set", fontsize=10)
plt.ylabel("Absolute Prediction Error", fontsize=10)
plt.legend(loc="upper right")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
ad_plot_path = os.path.join(
    output_dir, f"{target_id}_Applicability_Domain_Plot.png"
)
plt.savefig(ad_plot_path, dpi=300)
plt.show()
print(f"-> Saved Applicability Domain Plot to: {ad_plot_path}")
# ==========================================
# Cell 14: Phytochemical Virtual Screening (Clean Output)
# ==========================================
phyto_file_path = os.path.join(output_dir, "Argemone mexicana.xlsx")
if os.path.exists(phyto_file_path):
    phyto_df = pd.read_excel(phyto_file_path)
    
    # Strip hidden whitespaces from columns
    phyto_df.columns = phyto_df.columns.str.strip()
    phyto_fps_bit = []
    valid_compounds = []
    for _, row in phyto_df.iterrows():
        smi_val = row["SMILES"] if "SMILES" in phyto_df.columns else row.get("Smile")
        if pd.isna(smi_val):
            continue
        mol = Chem.MolFromSmiles(str(smi_val))
        if mol is not None:
            phyto_fps_bit.append(np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)))
            valid_compounds.append(row["Phytochemical name"])
    if len(valid_compounds) > 0:
        phyto_X = np.array(phyto_fps_bit)
        phyto_pred_pIC50 = rf_model.predict(phyto_X)
        phyto_predictions = pd.DataFrame({
            "Compound": valid_compounds,
            "Predicted_pIC50": phyto_pred_pIC50
        }).sort_values(by="Predicted_pIC50", ascending=False).reset_index(drop=True)
        print(f"\n--- Sorted Phytochemical {target_id} Predictions ---")
        print(phyto_predictions.to_string(index=False))
        pred_output_path = os.path.join(output_dir, f"{target_id}_Sorted_Phytochemical_Predictions.xlsx")
        phyto_predictions.to_excel(pred_output_path, index=False)
        print(f"\n-> Successfully saved screening sheet to: {pred_output_path}")
    else:
        print("No valid chemical structures found in the phytochemical file.")
else:
    print(f"Phytochemical file not found at: {phyto_file_path}")
