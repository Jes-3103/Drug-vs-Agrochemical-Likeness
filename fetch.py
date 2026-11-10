# fetch.py
# Fetches SMILES from PubChem for each compound name in data/compounds.csv
# and saves them to data/compounds_with_smiles.csv. Needs an internet connection.

import time
import pandas as pd
import pubchempy as pcp

df = pd.read_csv("data/compounds.csv")

smiles_list = []
for name in df["name"]:
    try:
        hits = pcp.get_compounds(name, "name")
        smiles = hits[0].smiles if hits else None
    except Exception as e:
        print(f"Could not fetch {name}: {e}")
        smiles = None
    print(f"{name}: {smiles}")
    smiles_list.append(smiles)
    time.sleep(0.3)   # be polite to PubChem's servers

df["smiles"] = smiles_list
df.to_csv("data/compounds_with_smiles.csv", index=False)
print("saved data/compounds_with_smiles.csv")
