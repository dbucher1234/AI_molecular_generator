#!/usr/bin/env python
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, AllChem
from rdkit.DataStructs import TanimotoSimilarity

# 1) Load your sampled SMILES + NLL
df = pd.read_csv("step3_samples.csv")

# 2) Reference diphenhydramine
dip_smiles = "CN(C)CCOC(c1ccccc1)c2ccccc2"
dip_mol    = Chem.MolFromSmiles(dip_smiles)
dip_fp     = AllChem.GetMorganFingerprintAsBitVect(dip_mol, radius=3, nBits=2048)
dip_logp   = Crippen.MolLogP(dip_mol)

# 3) Compute scores
results = []
for _, row in df.iterrows():
    mol = Chem.MolFromSmiles(row["SMILES"])
    if not mol:
        continue

    logp = Crippen.MolLogP(mol)
    fp   = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048)
    sim  = TanimotoSimilarity(dip_fp, fp)

    # continuous reward for lower logP
    logp_score = max(0.0, min(1.0, (dip_logp - logp) / dip_logp)) if dip_logp else 0.0
    combined   = 0.5 * sim + 0.5 * logp_score

    results.append({
        "SMILES":       row["SMILES"],
        "NLL":          row["NLL"],
        "logP":         round(logp, 2),
        "Similarity":   round(sim, 3),
        "logP_score":   round(logp_score, 3),
        "Combined":     round(combined, 3)
    })

# 4) Build DataFrame, sort, then apply hard filters
df_out = pd.DataFrame(results)
df_out = df_out.sort_values("Combined", ascending=False)

# Hard‐filters:
df_filtered = df_out[
    (df_out["Similarity"] >= 0.3) &
    (df_out["logP"] <= dip_logp)
]

# 5) Save and print top hits
df_filtered.to_csv("step4_rdkit_filtered.csv", index=False)
print(df_filtered.head(10))

