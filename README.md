# ProCMT-QSAR — Protein-Conditioned Multimodal Multi-Task QSAR for Endocrine Toxicity

Reproducibility repository for:

> **ProCMT-QSAR: An Interpretable Target-Aware Protein-Conditioned Multimodal Multi-Task QSAR Model for Toxicity Prediction of Endocrine-Disrupting Chemicals**
>
> Fauzan Syarif Nursyafi<sup>1</sup>, Yunendah Nur Fuadah<sup>2</sup>, Ulfa Latifa Hanum<sup>3</sup>, Ki Moo Lim<sup>1,4,5,*</sup>
>
> <sup>1</sup>Department of Biomedical Engineering, Kumoh National Institute of Technology, Gumi, 39177, Republic of Korea
> <sup>2</sup>Telecommunication Engineering Study Program, School of Electrical Engineering, Telkom University, Jl. Telekomunikasi No. 1, Bandung 40257, West Java, Indonesia
> <sup>3</sup>Department of IT Convergence Engineering, Kumoh National Institute of Technology, Gumi, 39177, Republic of Korea
> <sup>4</sup>Meta Heart Co., Ltd, Gumi, 39253, Republic of Korea
> <sup>5</sup>High-tech Medical Equipment Research Institute, Kumoh National Institute of Technology, Gumi, 39177, Republic of Korea
>
> \*Corresponding author: Ki Moo Lim, [kmlim@kumoh.ac.kr](mailto:kmlim@kumoh.ac.kr)

---

## Overview

This repo provides a step-by-step tutorial pipeline that reproduces the complete methodology from raw ChEMBL data to model predictions and explainability analysis.

The model predicts whether a compound is **active** (IC50 ≤ 1,000 nM) or **non-active** (IC50 ≥ 10,000 nM) against six endocrine-relevant protein targets simultaneously, using a protein-conditioned multimodal multi-task architecture.

| Target | Biological role |
|--------|----------------|
| AR (Androgen receptor) | Androgen signaling |
| CYP19A1 (Aromatase) | Estrogen biosynthesis |
| ESR1-α (Estrogen receptor α) | Estrogen receptor signaling |
| ESR2-β (Estrogen receptor β) | Estrogen receptor signaling |
| HSD11B1 (11β-HSD1) | Glucocorticoid metabolism |
| PPARγ (PPARG) | Metabolic nuclear receptor regulation |

---

## ProCMT-QSAR Framework

![ProCMT-QSAR Framework](ProCMTQSAR_framework.png)

The framework integrates six complementary molecular representations (fingerprints, physicochemical descriptors, Mordred 2D descriptors, molecular graph, SMILES sequence embedding, and 3D conformer descriptors) with protein-conditioned cross-attention and a modality-gated Mixture-of-Experts fusion module, trained jointly across six endocrine-related endpoints under a multi-task learning scheme.

---

## Repository structure

```
ProCMT-QSAR/
├── README.md
├── LICENSE
├── requirements.txt
├── config.yaml                        ← central configuration (paths, thresholds, seeds)
├── ProCMTQSAR_framework.png
│
├── data/
│   ├── raw/README.md                  ← ChEMBL download instructions
│   ├── processed/                     ← curated data output
│   └── demo/                          ← small demo subset for quick-run
│
├── utils/
│   ├── __init__.py
│   ├── chem_utils.py                  ← shared chemistry functions
│   └── split_utils.py                 ← scaffold split and leakage check
│
└── notebooks/
    ├── 01_dataset_preprocessing.ipynb
    ├── 02_chemical_space_and_ad.ipynb
    ├── 03_molecular_representation.ipynb
    ├── 04_protein_embedding.ipynb
    ├── 05_model_training.ipynb
    ├── 06_performance_evaluation.ipynb
    └── 07_explainable_ai.ipynb
```

> **Dataset & pretrained checkpoints:** the curated dataset and pretrained model checkpoints used in this study can be downloaded from our Zenodo repository (DOI: [10.5281/zenodo.dummy](https://doi.org/10.5281/zenodo.dummy)).

---

## Quick start

### 1. Clone the repo

```bash
git clone https://github.com/kit-cml/ProCMT-QSAR.git
cd ProCMT-QSAR
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

Follow the instructions in `data/raw/README.md` to download the six ChEMBL IC50 CSV files.

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
| 03 Molecular representation | Six modalities: FP, descriptors, graph, SMILES, 3D | Optional | ~30–60 min |
| 04 Protein embedding | ESM-2 embeddings for six targets | Optional | ~5 min |
| 05 Model training | Multimodal multi-task QSAR, 10-fold CV | Yes | ~12–24 hr |
| 06 Performance evaluation | Bootstrap CI, baselines, ROC/PR curves | Optional | ~30 min |
| 07 Explainable AI | SHAP, fragment, atom attribution, gate weights | Optional | ~1–2 hr |

---

## Key methodological decisions

- **Labeling:** Clean-margin binary labels — active ≤ 1,000 nM, non-active ≥ 10,000 nM, intermediate removed
- **Split:** Stratified scaffold split (Bemis–Murcko) — 80% train pool, 20% held-out test
- **Evaluation:** 10-fold stratified scaffold CV within train pool; held-out test used exactly once
- **Balancing:** Random undersampling of majority class, max 3:1 ratio (train pool only)
- **Architecture:** Protein-conditioned cross-attention + modality-gated Mixture-of-Experts fusion
- **Metrics:** AUROC, AUPRC, MCC, Balanced Accuracy, Sensitivity, Specificity, F1

---

## Citation

```bibtex
@article{nursyafi2026procmtqsar,
  title   = {ProCMT-QSAR: An Interpretable Target-Aware Protein-Conditioned
             Multimodal Multi-Task QSAR Model for Toxicity Prediction of
             Endocrine-Disrupting Chemicals},
  author  = {Nursyafi, Fauzan Syarif and Fuadah, Yunendah Nur and
             Hanum, Ulfa Latifa and Lim, Ki Moo},
  journal = {Journal Name},
  year    = {2026},
}
```

---

## License

MIT License — see `LICENSE` for details.
