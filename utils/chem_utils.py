"""
utils/chem_utils.py
-------------------
Shared chemistry utility functions used across all notebooks in the
QSAR Endocrine Toxicity pipeline.

Functions
---------
get_canonical_smiles      : Parse and canonicalize a SMILES string
get_largest_fragment      : Keep only the largest organic fragment (for salts)
get_bemis_murcko_scaffold : Compute Bemis-Murcko scaffold (SMILES string)
convert_to_nM             : Convert IC50 value to nanomolar units
compute_pIC50             : Compute pIC50 = 9 - log10(IC50_nM)
normalize_relation        : Standardize Standard Relation strings from ChEMBL
assign_binary_label       : Apply clean-margin binary labeling (IC50-based)
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger

# Suppress RDKit warnings globally
RDLogger.DisableLog("rdApp.*")


# ── Unit conversion map (to nM) ───────────────────────────────────────────────
UNIT_TO_NM = {
    "nM": 1.0,
    "nm": 1.0,
    "uM": 1_000.0,
    "um": 1_000.0,
    "µM": 1_000.0,
    "µm": 1_000.0,
    "mM": 1_000_000.0,
    "mm": 1_000_000.0,
    "M":  1_000_000_000.0,
    "pM": 0.001,
    "pm": 0.001,
}


def get_canonical_smiles(smiles: str) -> str | None:
    """
    Parse a SMILES string and return its canonical form.
    Returns None if the SMILES is invalid.

    Parameters
    ----------
    smiles : str
        Input SMILES string.

    Returns
    -------
    str or None
        Canonical SMILES, or None if invalid.
    """
    if not isinstance(smiles, str) or smiles.strip() == "":
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def get_largest_fragment(smiles: str) -> str | None:
    """
    For salts and multi-component entries, retain only the largest
    organic fragment. Returns the canonical SMILES of that fragment.

    Parameters
    ----------
    smiles : str
        Input SMILES string (may contain dots for mixtures/salts).

    Returns
    -------
    str or None
        Canonical SMILES of the largest fragment, or None if invalid.
    """
    if not isinstance(smiles, str) or smiles.strip() == "":
        return None

    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None

    # Split on fragments and pick the largest by heavy atom count
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None

    largest = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    return Chem.MolToSmiles(largest, canonical=True)


def get_bemis_murcko_scaffold(smiles: str) -> str | None:
    """
    Compute the Bemis-Murcko scaffold for a given canonical SMILES.

    Acyclic molecules produce an empty scaffold string — these are assigned
    the canonical SMILES itself as a fallback key so they don't all collapse
    into a single 'empty scaffold' group.

    Parameters
    ----------
    smiles : str
        Canonical SMILES string.

    Returns
    -------
    str or None
        Scaffold SMILES string, or None if parsing fails.
    """
    if not isinstance(smiles, str):
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    scaffold = MurckoScaffold.MurckoScaffoldSmiles(
        mol=mol, includeChirality=False
    )

    # Fallback for acyclic molecules (empty scaffold string)
    if scaffold == "" or scaffold is None:
        return smiles  # use molecule itself as unique scaffold key

    return scaffold


def convert_to_nM(value: float, unit: str) -> float | None:
    """
    Convert an IC50 value from its original unit to nanomolar (nM).

    Parameters
    ----------
    value : float
        Numeric IC50 value.
    unit : str
        Unit string from ChEMBL (e.g. 'nM', 'uM', 'µM').

    Returns
    -------
    float or None
        IC50 in nM, or None if unit is not recognized.
    """
    if not isinstance(value, (int, float)) or np.isnan(value):
        return None
    if value <= 0:
        return None

    multiplier = UNIT_TO_NM.get(str(unit).strip())
    if multiplier is None:
        return None

    return float(value) * multiplier


def compute_pIC50(ic50_nM: float) -> float | None:
    """
    Compute pIC50 = 9 - log10(IC50_nM).

    This transforms IC50 in nM to the negative log10 scale relative to 1 M:
        pIC50 = -log10(IC50_M) = -log10(IC50_nM * 1e-9) = 9 - log10(IC50_nM)

    Parameters
    ----------
    ic50_nM : float
        IC50 value in nanomolar (nM). Must be positive.

    Returns
    -------
    float or None
        pIC50 value, or None if ic50_nM is invalid.
    """
    if ic50_nM is None or not isinstance(ic50_nM, (int, float)):
        return None
    if np.isnan(ic50_nM) or ic50_nM <= 0:
        return None
    return 9.0 - np.log10(float(ic50_nM))


def normalize_relation(relation: str) -> str:
    """
    Standardize the Standard Relation field from ChEMBL.

    ChEMBL may export relations as '=', '>', '<', '>=', '<=', or with
    surrounding whitespace/quotes. This returns a clean string.

    Parameters
    ----------
    relation : str
        Raw Standard Relation value from ChEMBL.

    Returns
    -------
    str
        Cleaned relation string, e.g. '=', '>', '<', '>=', '<=', or 'unknown'.
    """
    if not isinstance(relation, str):
        return "unknown"

    r = relation.strip().strip('"').strip("'")

    mapping = {
        "=":  "=",
        ">":  ">",
        "<":  "<",
        ">=": ">=",
        "<=": "<=",
        "~":  "~",
    }
    return mapping.get(r, "unknown")


def assign_binary_label(
    ic50_nM: float,
    active_threshold: float = 1000.0,
    inactive_threshold: float = 10000.0,
) -> int | None:
    """
    Assign a binary label using the clean-margin IC50 thresholding scheme.

    Label rules:
        active   (1) : IC50 <= active_threshold    (default: 1000 nM)
        non-active(0) : IC50 >= inactive_threshold  (default: 10000 nM)
        removed (None): intermediate zone (1000 < IC50 < 10000 nM)

    Parameters
    ----------
    ic50_nM : float
        IC50 value in nanomolar.
    active_threshold : float
        Upper bound for active label (inclusive). Default: 1000 nM.
    inactive_threshold : float
        Lower bound for non-active label (inclusive). Default: 10000 nM.

    Returns
    -------
    int or None
        1 (active), 0 (non-active), or None (intermediate — to be removed).
    """
    if ic50_nM is None or not isinstance(ic50_nM, (int, float)):
        return None
    if np.isnan(ic50_nM) or ic50_nM <= 0:
        return None

    if ic50_nM <= active_threshold:
        return 1
    if ic50_nM >= inactive_threshold:
        return 0
    return None  # intermediate zone — will be removed


def assign_removed_reason(row) -> str:
    """
    Return a human-readable reason for why a record was removed.
    Intended to be used with DataFrame.apply() on the removed subset.

    Parameters
    ----------
    row : pd.Series
        A row from a preprocessed DataFrame containing at minimum
        'IC50_nM' and 'Standard Relation Clean'.

    Returns
    -------
    str
        Short description of the removal reason.
    """
    ic50 = row.get("IC50_nM", None)
    relation = row.get("Standard Relation Clean", None)

    if ic50 is None or (isinstance(ic50, float) and np.isnan(ic50)):
        return "missing_or_invalid_IC50"
    if relation not in ("=", None, "unknown"):
        return f"non_exact_relation ({relation})"
    if 1000 < ic50 < 10000:
        return "intermediate_zone (1000 < IC50 < 10000 nM)"
    return "other"
