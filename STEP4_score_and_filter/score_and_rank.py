#!/usr/bin/env python
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, AllChem
from rdkit.DataStructs import TanimotoSimilarity

# 1) Load sampled SMILES + NLL
df = pd.read_csv("step3_samples.csv")

# 2) Reference diphenhydramine
dip_smiles = "CN(C)CCOC(c1ccccc1)c2ccccc2"
dip_mol = Chem.MolFromSmiles(dip_smiles)
dip_fp = AllChem.GetMorganFingerprintAsBitVect(dip_mol, radius=3, nBits=2048)
dip_logp = Crippen.MolLogP(dip_mol)

# 3) Determine NLL normalization range
nll_min, nll_max = df["NLL"].min(), df["NLL"].max()

# 4) Compute filtered & scored list
results = []
for _, row in df.iterrows():
    smi = row["SMILES"]
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        continue

    # Hard filter on logP
    logp = Crippen.MolLogP(mol)
    if logp > dip_logp:
        continue

    # Compute similarity
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048)
    similarity = TanimotoSimilarity(dip_fp, fp)

    # Reverse-normalized NLL
    nll = row["NLL"]
    reverse_nll = 1 - ((nll - nll_min) / (nll_max - nll_min)) if nll_max > nll_min else 1.0

    # Combined similarity score
    combined = (reverse_nll + similarity) / 2

    results.append({
        "SMILES": smi,
        "logP": round(logp, 2),
        "Similarity": round(similarity, 3),
        "Reverse_NLL": round(reverse_nll, 3),
        "Combined": round(combined, 3)
    })

# 5) Build DataFrame and sort
df_final = pd.DataFrame(results).sort_values("Combined", ascending=False)

# 6) Save output
df_final.to_csv("step4_final.csv", index=False)
print(df_final.head(10))

