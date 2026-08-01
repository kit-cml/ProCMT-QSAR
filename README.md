# Protein-Conditioned Multimodal Multi-Task QSAR — Endocrine Toxicity

Reproducibility repository for:

> **Protein-Aware Multimodal Multi-Task QSAR Classification of IC50-Based Endocrine Toxicity Mechanism Endpoints across AR, CYP19A1, ESR1, ESR2, HSD11B1, and PPARG under Stratified Scaffold Evaluation**

---

## Overview

This repo provides a step-by-step tutorial pipeline that reproduces the complete methodology from raw ChEMBL data to model predictions and explainability analysis.

The model predicts whether a compound is **active** (IC50 ≤ 1,000 nM) or **non-active** (IC50 ≥ 10,000 nM) against six endocrine-relevant protein targets simultaneously, using a protein-conditioned multimodal multi-task architecture.

| Target | Biological role |
|--------|----------------|
| AR (Androgen receptor) | Androgen signaling |
| CYP19A1 (Aromatase) | Estrogen biosynthesis |
| ESR1-alpha (Estrogen receptor α) | Estrogen receptor signaling |
| ESR2-beta (Estrogen receptor β) | Estrogen receptor signaling |
| HSD11B1 (11β-HSD1) | Glucocorticoid metabolism |
| PPARG (PPARγ) | Metabolic nuclear receptor regulation |

---

## Repository structure

```
qsar-endocrine-toxicity/
├── README.md
├── LICENSE
├── requirements.txt
├── config.yaml                        ← central configuration (paths, thresholds, seeds)
│
├── data/
│   ├── raw/README.md                  ← ChEMBL download instructions
│   ├── processed/                     ← curated data output
│   └── demo/                          ← small demo subset for quick-run
│
├── utils/
│   ├── chem_utils.py                  ← shared chemistry functions
│   ├── split_utils.py                 ← scaffold split and leakage check
│   └── eval_utils.py                  ← metrics and bootstrap CI
│
├── notebooks/
│   ├── 01_dataset_preprocessing.ipynb
│   ├── 02_chemical_space_and_ad.ipynb
│   ├── 03_molecular_representation.ipynb
│   ├── 04_protein_embedding.ipynb
│   ├── 05_model_training.ipynb
│   ├── 06_performance_evaluation.ipynb
│   └── 07_explainable_ai.ipynb
│
├── scripts/                           ← full-scale GPU scripts (see below)
│   ├── run_full_training.py
│   └── run_full_evaluation.py
│
└── checkpoints/README.md             ← link to pretrained model weights
```

---

## Quick start

### 1. Clone the repo

```bash
git clone https://github.com/your-username/qsar-endocrine-toxicity.git
cd qsar-endocrine-toxicity
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

For GPU support (recommended for notebooks 05–07), install PyTorch with CUDA from [pytorch.org](https://pytorch.org/get-started/locally/) first, then run the above.

### 3. Configure paths

Open `config.yaml` and set `raw_data_root` and `processed_root` to match your local machine:

```yaml
paths:
  raw_data_root:  "path/to/your/raw/chembl/files"
  processed_root: "path/to/your/output/folder"
```

### 4. Download raw data

Follow the instructions in `data/raw/README.md` to download the 6 ChEMBL IC50 CSV files.

### 5. Run the pipeline

Run notebooks in order:

```
01 → 02 → 03 → 04 → 05 → 06 → 07
```

Each notebook is self-contained and loads its inputs from the outputs of the previous one.

---

## Demo mode (no data download needed)

Set `demo.enabled: true` in `config.yaml` to run the entire pipeline on a small synthetic dataset (~200 molecules per target). This takes about 5 minutes on CPU and requires no GPU.

---

## Notebook overview

| Notebook | Description | GPU? | Time (full) |
|----------|-------------|------|-------------|
| 01 Dataset preprocessing | Curation, labeling, scaffold split | No | ~5–10 min |
| 02 Chemical space & AD | UMAP, Tanimoto, scaffold diversity | No | ~5 min |
| 03 Molecular representation | 6 modalities: FP, descriptors, graph, SMILES, 3D | Optional | ~30–60 min |
| 04 Protein embedding | ESM-2 embeddings for 6 targets | Optional | ~5 min |
| 05 Model training | Multimodal multi-task QSAR, 10-fold CV | Yes | ~12–24 hr |
| 06 Performance evaluation | Bootstrap CI, baselines, ROC/PR curves | Optional | ~30 min |
| 07 Explainable AI | SHAP, fragment, atom attribution, gate weights | Optional | ~1–2 hr |

---

## Key methodological decisions

- **Labeling:** Clean-margin binary labels — active ≤ 1,000 nM, non-active ≥ 10,000 nM, intermediate removed
- **Split:** Stratified scaffold split (Bemis-Murcko) — 80% train pool, 20% held-out test
- **Evaluation:** 10-fold stratified scaffold CV within train pool; held-out test used exactly once
- **Balancing:** Random undersampling of majority class, max 3:1 ratio (train pool only)
- **Architecture:** Protein-conditioned cross-attention + modality-gated Mixture-of-Experts fusion
- **Metrics:** AUROC, AUPRC, MCC, Balanced Accuracy, Sensitivity, Specificity, F1

---

## Citation

```bibtex
@article{yourname2026qsar,
  title   = {Protein-Aware Multimodal Multi-Task QSAR Classification of IC50-Based
             Endocrine Toxicity Mechanism Endpoints across AR, CYP19A1, ESR1,
             ESR2, HSD11B1, and PPARG under Stratified Scaffold Evaluation},
  author  = {Author, A. and Author, B.},
  journal = {Journal Name},
  year    = {2026},
}
```

---

## License

MIT License — see `LICENSE` for details.
