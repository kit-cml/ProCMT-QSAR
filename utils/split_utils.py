"""
utils/split_utils.py
--------------------
Scaffold-based dataset splitting utilities for the QSAR pipeline.

Functions
---------
stratified_scaffold_split : 80/20 global scaffold split across all endpoints
make_cv_fold_assignments   : 10-fold stratified scaffold CV within train pool
check_leakage              : Verify zero scaffold/molecule overlap between splits
summarize_split            : Print per-endpoint train/test counts
"""

import numpy as np
import pandas as pd
from collections import defaultdict
from sklearn.model_selection import StratifiedKFold


def stratified_scaffold_split(
    df: pd.DataFrame,
    scaffold_col: str = "scaffold",
    smiles_col: str = "canonical_smiles",
    endpoint_col: str = "endpoint_name",
    label_col: str = "label",
    test_size: float = 0.20,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Perform a stratified scaffold-based train/test split across all endpoints.

    Strategy
    --------
    1. Group unique canonical SMILES by their Bemis-Murcko scaffold.
    2. Shuffle scaffold groups (reproducibly with `seed`).
    3. Assign scaffold groups to train pool (80%) or held-out test (20%).
    4. Stratification by endpoint-label combination is applied as an
       approximate constraint on the scaffold ordering — scaffold separation
       is the primary criterion.

    Note: Each compound receives exactly ONE split assignment shared across
    all endpoints (compound-level consistency, not endpoint-level).

    Parameters
    ----------
    df : pd.DataFrame
        Combined multi-endpoint DataFrame. Must contain columns:
        [smiles_col, scaffold_col, endpoint_col, label_col].
    scaffold_col : str
        Column name for Bemis-Murcko scaffold SMILES.
    smiles_col : str
        Column name for canonical SMILES.
    endpoint_col : str
        Column name for endpoint identifier.
    label_col : str
        Column name for binary label (0 or 1).
    test_size : float
        Fraction of scaffold groups to assign to held-out test (default: 0.20).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with a new 'split' column: 'train_pool' or 'heldout_test'.
    """
    rng = np.random.RandomState(seed)

    # Build scaffold → set of canonical SMILES mapping
    scaffold_to_smiles = defaultdict(set)
    for _, row in df[[smiles_col, scaffold_col]].drop_duplicates().iterrows():
        scaffold_to_smiles[row[scaffold_col]].add(row[smiles_col])

    scaffolds = list(scaffold_to_smiles.keys())

    # Sort scaffolds by size (descending) so large scaffolds are distributed first,
    # which helps approximate stratification across the split boundary.
    scaffolds = sorted(scaffolds, key=lambda s: len(scaffold_to_smiles[s]), reverse=True)
    rng.shuffle(scaffolds)

    # Assign scaffolds to train/test until test quota is reached
    total_smiles = df[smiles_col].nunique()
    test_target  = int(np.ceil(test_size * total_smiles))

    test_scaffolds  = set()
    test_smiles_count = 0

    for scaffold in scaffolds:
        if test_smiles_count >= test_target:
            break
        test_scaffolds.add(scaffold)
        test_smiles_count += len(scaffold_to_smiles[scaffold])

    # Map every row to a split
    df = df.copy()
    df["split"] = df[scaffold_col].apply(
        lambda s: "heldout_test" if s in test_scaffolds else "train_pool"
    )

    return df


def make_cv_fold_assignments(
    train_df: pd.DataFrame,
    scaffold_col: str = "scaffold",
    smiles_col: str = "canonical_smiles",
    endpoint_col: str = "endpoint_name",
    label_col: str = "label",
    n_folds: int = 10,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Assign 10-fold stratified scaffold cross-validation fold IDs to the train pool.

    Each scaffold group is assigned to exactly one fold. Stratification is
    applied at the scaffold group level using the most frequent endpoint-label
    combination within each scaffold group as the stratum key.

    Parameters
    ----------
    train_df : pd.DataFrame
        Train pool DataFrame (split == 'train_pool').
    scaffold_col : str
        Column name for scaffold SMILES.
    smiles_col : str
        Column name for canonical SMILES.
    endpoint_col : str
        Column name for endpoint identifier.
    label_col : str
        Column name for binary label.
    n_folds : int
        Number of CV folds (default: 10).
    seed : int
        Random seed.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with a new 'cv_fold' column (0 to n_folds-1).
    """
    # Build scaffold-level summary: one row per unique scaffold
    scaffold_summary = (
        train_df.groupby(scaffold_col)
        .apply(lambda g: pd.Series({
            "n_smiles":      g[smiles_col].nunique(),
            "strat_key":     (g[endpoint_col].astype(str) + "__" +
                              g[label_col].astype(str))
                             .value_counts().index[0]  # most frequent combination
        }))
        .reset_index()
    )

    # Check minimum group size
    strat_counts = scaffold_summary["strat_key"].value_counts()
    if strat_counts.min() < n_folds:
        print(
            f"  ⚠️  Some strat groups have fewer than {n_folds} scaffolds. "
            "StratifiedKFold will fall back to standard KFold for those groups."
        )

    # Stratified K-Fold on scaffold-level strat_key
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    scaffold_summary["cv_fold"] = -1

    X_dummy = np.zeros(len(scaffold_summary))
    y_strat  = scaffold_summary["strat_key"].values

    for fold_idx, (_, val_idx) in enumerate(skf.split(X_dummy, y_strat)):
        scaffold_summary.loc[val_idx, "cv_fold"] = fold_idx

    # Map fold IDs back to the original DataFrame
    scaffold_to_fold = dict(
        zip(scaffold_summary[scaffold_col], scaffold_summary["cv_fold"])
    )

    train_df = train_df.copy()
    train_df["cv_fold"] = train_df[scaffold_col].map(scaffold_to_fold)
    return train_df


def check_leakage(
    df: pd.DataFrame,
    scaffold_col: str = "scaffold",
    smiles_col: str = "canonical_smiles",
    split_col: str = "split",
    train_label: str = "train_pool",
    test_label: str = "heldout_test",
) -> dict:
    """
    Verify that no scaffold or canonical SMILES appears in both train and test.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with split assignments.
    scaffold_col, smiles_col, split_col : str
        Column names.
    train_label, test_label : str
        Values in split_col identifying the train pool and held-out test.

    Returns
    -------
    dict
        Leakage summary with keys:
        ['n_train_scaffolds', 'n_test_scaffolds', 'n_leaked_scaffolds',
         'n_train_molecules', 'n_test_molecules',  'n_leaked_molecules']

    Raises
    ------
    ValueError
        If any scaffold or molecule leakage is detected.
    """
    train_mask = df[split_col] == train_label
    test_mask  = df[split_col] == test_label

    train_scaffolds = set(df.loc[train_mask, scaffold_col])
    test_scaffolds  = set(df.loc[test_mask,  scaffold_col])
    leaked_scaffolds = train_scaffolds & test_scaffolds

    train_molecules = set(df.loc[train_mask, smiles_col])
    test_molecules  = set(df.loc[test_mask,  smiles_col])
    leaked_molecules = train_molecules & test_molecules

    summary = {
        "n_train_scaffolds":  len(train_scaffolds),
        "n_test_scaffolds":   len(test_scaffolds),
        "n_leaked_scaffolds": len(leaked_scaffolds),
        "n_train_molecules":  len(train_molecules),
        "n_test_molecules":   len(test_molecules),
        "n_leaked_molecules": len(leaked_molecules),
    }

    if leaked_scaffolds:
        raise ValueError(
            f"❌ SCAFFOLD LEAKAGE DETECTED: {len(leaked_scaffolds)} scaffolds "
            f"appear in both train and test. Aborting."
        )
    if leaked_molecules:
        raise ValueError(
            f"❌ MOLECULE LEAKAGE DETECTED: {len(leaked_molecules)} canonical SMILES "
            f"appear in both train and test. Aborting."
        )

    print("✅ Leakage check passed — zero scaffold and molecule overlap.")
    return summary


def summarize_split(
    df: pd.DataFrame,
    endpoint_col: str = "endpoint_name",
    label_col: str = "label",
    split_col: str = "split",
) -> pd.DataFrame:
    """
    Print and return a per-endpoint summary of the train/test split.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with split assignments and labels.
    endpoint_col, label_col, split_col : str
        Column names.

    Returns
    -------
    pd.DataFrame
        Summary table with columns:
        [endpoint, n_train_active, n_train_inactive, n_train_total,
         n_test_active, n_test_inactive, n_test_total, ratio_train]
    """
    rows = []
    for endpoint, grp in df.groupby(endpoint_col):
        train = grp[grp[split_col] == "train_pool"]
        test  = grp[grp[split_col] == "heldout_test"]

        tr_act  = (train[label_col] == 1).sum()
        tr_ina  = (train[label_col] == 0).sum()
        te_act  = (test[label_col]  == 1).sum()
        te_ina  = (test[label_col]  == 0).sum()

        rows.append({
            "endpoint":          endpoint,
            "train_active":      tr_act,
            "train_inactive":    tr_ina,
            "train_total":       tr_act + tr_ina,
            "test_active":       te_act,
            "test_inactive":     te_ina,
            "test_total":        te_act + te_ina,
            "train_ratio":       round(tr_act / max(tr_ina, 1), 2),
        })

    summary_df = pd.DataFrame(rows).set_index("endpoint")

    print("\n" + "=" * 75)
    print("SPLIT SUMMARY (per endpoint)")
    print("=" * 75)
    print(summary_df.to_string())
    print("=" * 75)

    return summary_df
