# AI Molecular Generator
Generating Polar Analogues of Diphenhydramine with REINVENT4

<p align="center">
  <img src="images/Diphen.png" alt="diphenhydramine (DHM) structure" />
</p>

## 📖 Introduction

This tutorial demonstrates how to use **REINVENT4** to generate de-novo small-molecule analogues of the classic antihistamine **diphenhydramine** (SMILES: `CN(C)CCOC(c1ccccc1)c2ccccc2`). Diphenhydramine is quite lipophilic (log P ≈ 3.4–3.7), which, along with its pKₐ and high passive permeability, drives strong CNS uptake—total brain levels ~18× plasma in rats, with an unbound‐drug brain/plasma ratio of ∼4–7. 

Our goal is to bias the generator toward more polar compounds that are less likely to cross the blood–brain barrier, **reducing central-nervous-system side effects**.

---

## 🔧 Prerequisites

- **REINVENT4** installed (see [REINVENT4 installation guide](https://github.com/MolecularAI/Reinvent4))  
- Python 3.8+  
- RDKit  
- A working CUDA-enabled GPU (optional, but recommended)  

---

## 🗂️ Repository Structure

```
.
├── STEP1_pretrain_and_transfer_learning/
├── STEP2_prepare_reference/
├── STEP3_sampling/
├── STEP4_score_and_filter/
├── STEP5_remove_existing/

```

---

## ⚙️ Workflow Overview

1. **Pretrain & Transfer Learning**  
   Fine-tune a REINVENT4 model on an antihistamine-focused dataset to bias it toward relevant chemotypes.

2. **Create Closely Related Molecules**  
   Use diphenhydramine as a starting point to bias the network toward closely related molecules.

3. **Sampling**  
   Generate a large pool of candidate analogues.

4. **Scoring & Filtering**  
   Compute predicted log P, CNS-MPO score, and scaffold novelty to select polar, BBB-impermeable scaffolds.

5. **Remove Known Molecules**  
   Eliminate any molecules already present in ChEMBL or PubChem to focus on novel hits.

---

## 🚀 STEP 1: Pretrain & Transfer Learning

### 🔍 What’s Happening in Transfer Learning?

1. **Starting from a Broad Prior**  
   We begin with `reinvent.prior`, a model pretrained on millions of drug-like molecules (from public databases, e.g. ChEMBL, PubChem) (Loeffler et al., 2024). This “prior” knows general medicinal-chemistry grammar: how to stitch atoms into plausible, synthesizable small molecules.

2. **Focusing on Antihistamines**  
   Next, we fine-tune (transfer-learn) that prior using our custom `all_antihistamines.smi` dataset—~300 known H₁-blockers from ChEMBL. This biases the model toward the scaffolds, functional groups, and chemotypes characteristic of antihistamines.

3. **Balancing Novelty vs. Familiarity**  
   We don’t want to lose all the broad-chemistry knowledge, nor merely memorize known antihistamines.  
   - **Similarity-pair filtering** enforces that, during training, generated molecules stay within a Tanimoto window (e.g. 0.7–1.0) of our antihistamine set.  
   - This encourages the model to explore **close analogues** of antihistamines, while the underlying prior still “remembers” general drug-like rules.

4. **Saving a New Specialized Prior**  
   The result is `my_project.prior`—a model that speaks both “general drug-design” (from the original prior) and “antihistamine” (from our fine-tuning).  
   You can now use it to sample entirely new, related scaffolds that retain key H₁-blocker motifs but may improve properties (e.g. polarity, hERG liability).

**Run:**  

   ```bash
   cd STEP1_pretrain_and_transfer_learning
   reinvent transfer_learning.toml
   ```

---

## Step 2: Mol2Mol Transfer Learning

In Step 2, we move from simply biasing our SMILES generator toward the general anti-histamines dataset (as in Step 1) to training a conditional model that learns to take an input scaffold and produce its close analogs. By fine-tuning on pairs of highly similar molecules, the Mol2Mol prior becomes specialized for lead optimization and analog design rather than broad, unconstrained sampling.

### Key Elements

- **Conditional Prior**  
  Start from the generic Mol2Mol model (`priors/mol2mol_scaffold_generic.prior`), which is architected to learn mappings between molecules.

- **Paired Training**  
  Supply pairs of SMILES (scaffold → analog) chosen by Tanimoto similarity (e.g. 0.99–1.0). The model learns to “translate” each scaffold into its near‐neighbor.

- **Focused Sampling**  
  After training, seeding the model with any scaffold yields close analogs, ideal for exploring small modifications around your lead compounds.


---

## 🎲 STEP 3: Sampling

Run the generator to produce a pool of candidates:

\`\`\`bash
reinvent4 sample \
  --model STEP1_pretrain_and_transfer_learning/final.ckpt \
  --prompts STEP2_prepare_reference/prompts.smi \
  --batch-size 128 \
  --samples 10000 \
  --output STEP3_sampling/raw_samples.smi
\`\`\`

---

## 📊 STEP 4: Scoring & Filtering

1. **Create Scoring Script**  
   Implement \`score_logp_cns.py\` to compute:
   - **Predicted log P** (target ≤ 2.5)  
   - **CNS-MPO** score (target ≥ 3)  
   - Tanimoto similarity to diphenhydramine

2. **Filter**  

   \`\`\`bash
   reinvent4 score \
     --input STEP3_sampling/raw_samples.smi \
     --scoring-func score_logp_cns.py \
     --thresholds '{"max_logp":2.5,"min_cns_mpo":3,"max_sim":0.8}' \
     --output STEP4_score_and_filter/filtered.smi
   \`\`\`

---

## 🚫 STEP 5: Remove Known Molecules

Remove any compounds already in public databases:

\`\`\`bash
reinvent4 dedupe \
  --input STEP4_score_and_filter/filtered.smi \
  --databases ChEMBL PubChem \
  --output STEP5_remove_existing/novel_hits.smi
\`\`\`

---

## 📈 Results & Next Steps

- Load \`novel_hits.smi\` into **RDKit**, **Molecular Notebook**, or **MOE**.  
- Cluster by Bemis–Murcko scaffold; pick diverse, polar scaffolds.  
- Prioritize top candidates for synthesis and in vitro testing.  
- Extend scoring to include predicted hERG liability, solubility, or synthetic accessibility.

---

## 🤝 Contributing

We welcome contributions! To help improve this workflow:

1. Fork the repository  
2. Add or refine scoring metrics (e.g., synthetic accessibility)  
3. Share your generated hits and any assay data  
4. Submit issues or pull requests  

---

## 📜 References

- O’Boyle, N.M. et al., *Mol2Mol perturbation methods*  
- Loeffler, H.H., He, J., Tibo, A. et al. Reinvent 4: Modern AI–driven generative molecule design. J Cheminform 16, 20 (2024). https://doi.org/10.1186/s13321-024-00812-5
