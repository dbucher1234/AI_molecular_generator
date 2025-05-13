from rdkit import Chem
from rdkit.Chem import Crippen, AllChem
from rdkit.DataStructs import TanimotoSimilarity
import pandas as pd

dip_smiles = "CN(C)CCOC(c1ccccc1)c2ccccc2"
dip_mol = Chem.MolFromSmiles(dip_smiles)
dip_fp = AllChem.GetMorganFingerprintAsBitVect(dip_mol, radius=3, nBits=2048)
dip_logp = Crippen.MolLogP(dip_mol)

df = pd.read_csv("step3_samples.csv")
results = []

for _, row in df.iterrows():
    mol = Chem.MolFromSmiles(row["SMILES"])
    if not mol:
        continue
    logp = Crippen.MolLogP(mol)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048)
    sim = TanimotoSimilarity(dip_fp, fp)
    logp_score = max(0.0, min(1.0, (dip_logp - logp) / dip_logp))
    combined = 0.5 * sim + 0.5 * logp_score
    results.append({
        "SMILES": row["SMILES"],
        "NLL": row["NLL"],
        "logP": round(logp, 2),
        "Similarity": round(sim, 3),
        "logP_score": round(logp_score, 3),
        "Combined": round(combined, 3)
    })

df_out = pd.DataFrame(results).sort_values("Combined", ascending=False)
df_out.to_csv("step4_rdkit_scored.csv", index=False)
print(df_out.head())

