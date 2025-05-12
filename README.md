# AI Molecular Generator
Generating Polar Analogues of Diphenhydramine with REINVENT4

<p align="center">
  <img src="images/Diphen.png" alt="diphenhydramine (DHM) structure" />
</p>

## 📖 Introduction

This tutorial demonstrates how to use **REINVENT4** to generate de-novo small-molecule analogues of the classic antihistamine **diphenhydramine** (SMILES: `CN(C)CCOC(c1ccccc1)c2ccccc2`).  
Our goal is to bias the generator toward more polar compounds that are less likely to cross the blood–brain barrier, reducing central-nervous-system side effects.

---

## 🔧 Prerequisites

- **REINVENT4** installed (see [REINVENT4 installation guide](https://github.com/MolecularAI/Reinvent4))  
- Python 3.8+  
- RDKit  
- A working CUDA-enabled GPU (optional, but highly recommended)  
- Basic familiarity with command-line operations  

\`\`\`bash
pip install reinvent4 rdkit-pypi
\`\`\`

---

## 🗂️ Repository Structure

\`\`\`
.
├── STEP1_pretrain_and_transfer_learning/
├── STEP2_prepare_reference/
├── STEP3_sampling/
├── STEP4_score_and_filter/
├── STEP5_remove_existing/
└── README.md
\`\`\`

---

## ⚙️ Workflow Overview

1. **Pretrain & Transfer Learning**  
   Fine-tune a REINVENT4 model on an antihistamine-focused dataset to bias it toward relevant chemotypes.

2. **Prepare Reference Pairs**  
   Convert diphenhydramine and simple perturbations (mol2mol variants) into tokenized prompts.

3. **Sampling**  
   Generate a large pool of candidate analogues.

4. **Scoring & Filtering**  
   Compute predicted log P, CNS-MPO score, and scaffold novelty to select polar, BBB-impermeable scaffolds.

5. **Remove Known Molecules**  
   Eliminate any molecules already present in ChEMBL or PubChem to focus on novel hits.

---

## 🚀 STEP 1: Pretrain & Transfer Learning

1. **Gather dataset**  
   - Collect ~1 000 known antihistamines (e.g., from ChEMBL).  
   - Format as one SMILES per line in \`data/antihistamines.smi\`.

2. **Configure transfer-learning**  
   Create \`step1_tl/config.yaml\`:

   \`\`\`yaml
   model:
     base_checkpoint: ./models/base_pretrained.ckpt
     target_smiles: ../data/antihistamines.smi
     epochs: 5
     learning_rate: 1e-4

   training:
     batch_size: 64
     device: cuda
   \`\`\`

3. **Run**  

   \`\`\`bash
   reinvent4 train \
     --config STEP1_pretrain_and_transfer_learning/step1_tl/config.yaml \
     --output-dir STEP1_pretrain_and_transfer_learning/
   \`\`\`

---

## 🛠️ STEP 2: Prepare Reference Pairs

1. **Generate Mol2Mol Variants**  
   - Swap phenyl → pyridyl, add small polar groups, etc.  
   - Save variants to \`STEP2_prepare_reference/mol2mol_variants.smi\`.

2. **Tokenize Prompts**  

   \`\`\`bash
   reinvent4 prepare \
     --input-smiles STEP2_prepare_reference/mol2mol_variants.smi \
     --output-dir STEP2_prepare_reference/ \
     --mode mol2mol
   \`\`\`

   Outputs \`prompts.smi\` ready for sampling.

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
- Zhavoronkov, A. et al., *REINVENT: de novo design* (J. Chem. Inf. Model., 2020)  
