# analyse.py
# Calculates descriptors and applies the drug- and agrochemical-likeness filters.
# Reads data/compounds_with_smiles.csv, writes results/screening_results.csv.

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

os.makedirs("results", exist_ok=True)


# --- Descriptors ---------------------------------------------------------

def calc_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return None
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),        # calculated (Crippen), not measured
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),   # counts N and O only
        "RotB": Lipinski.NumRotatableBonds(mol),
    }


df = pd.read_csv("data/compounds_with_smiles.csv")
desc = df["smiles"].apply(calc_descriptors)
print("Could not parse:", list(df[desc.isna()]["name"]))
df = df[desc.notna()].reset_index(drop=True)
df = pd.concat([df, pd.DataFrame(desc.dropna().tolist())], axis=1)


# --- Lipinski Rule of 5: at most one violation allowed --------------------

def lipinski_violations(row):
    return sum([row["MW"] > 500, row["LogP"] > 5,
                row["HBD"] > 5, row["HBA"] > 10])


df["Ro5_violations"] = df.apply(lipinski_violations, axis=1)
df["Lipinski_pass"] = df["Ro5_violations"] <= 1


# --- Veber rules for oral bioavailability ---------------------------------

df["Veber_pass"] = (df["RotB"] <= 10) & (df["TPSA"] <= 140)


# --- Tice rules for agrochemicals (Tice, Pest Manag. Sci. 2001) -----------
# These are ranges, not just upper limits, so zero violations is the pass mark.

def tice_violations(row, kind):
    if kind == "herbicide":
        checks = [not (150 <= row["MW"] <= 500), row["LogP"] > 3.5,
                  row["HBD"] > 3, not (2 <= row["HBA"] <= 12), row["RotB"] >= 12]
    else:
        checks = [not (150 <= row["MW"] <= 500), not (0 <= row["LogP"] <= 5),
                  row["HBD"] > 2, not (1 <= row["HBA"] <= 8), row["RotB"] >= 12]
    return sum(checks)


df["Tice_herb_viol"] = df.apply(tice_violations, axis=1, kind="herbicide")
df["Tice_ins_viol"] = df.apply(tice_violations, axis=1, kind="insecticide")
df["Tice_herb_pass"] = df["Tice_herb_viol"] == 0
df["Tice_ins_pass"] = df["Tice_ins_viol"] == 0


# --- PAINS: substructures that often give false positives in assays -------

params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
pains = FilterCatalog(params)


def pains_alert(smiles):
    entry = pains.GetFirstMatch(Chem.MolFromSmiles(smiles))
    return entry.GetDescription() if entry else ""


df["PAINS_alert"] = df["smiles"].apply(pains_alert)


# --- Save and summarise ---------------------------------------------------

df.round(2).to_csv("results/screening_results.csv", index=False)
print("saved results/screening_results.csv")

print()
print(df.groupby("category")[["MW", "LogP", "TPSA"]].mean().round(1))
print()
print(df.groupby("category")[["Lipinski_pass", "Tice_herb_pass",
                              "Tice_ins_pass"]].mean().round(2))
