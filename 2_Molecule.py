import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, Crippen, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from mordred import Calculator, descriptors as _mordred_desc

from utils.shared import (
    render_feature_table, offer_csv_download, offer_excel_download,
    pick_column, read_uploaded_table, save_fig_bytes,
)

st.set_page_config(page_title="Molecule Features", page_icon="⚛️", layout="wide")

plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 600, "savefig.bbox": "tight", "font.size": 11})

SMILES_COLS = ["smiles", "smile", "structure", "canonical_smiles", "smiles_string"]


@st.cache_resource
def get_mordred_calculator():
    return Calculator(_mordred_desc, ignore_3D=True)  # 1613 2-D descriptors


def analyze_smiles(s):
    mol = Chem.MolFromSmiles(str(s).strip())
    if mol is None:
        st.error("Could not parse that SMILES - check the string.")
        return
    mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
    hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    render_feature_table([
        ("Canonical SMILES", Chem.MolToSmiles(mol), "RDKit-normalised structure"),
        ("Formula", rdMolDescriptors.CalcMolFormula(mol), "Molecular formula"),
        ("Molecular weight", round(mw, 2), "g/mol"),
        ("LogP (Crippen)", round(logp, 2), "Lipophilicity"),
        ("TPSA", round(rdMolDescriptors.CalcTPSA(mol), 2), "Polar surface area (Å²)"),
        ("H-bond donors", hbd, "Lipinski rule: <= 5"),
        ("H-bond acceptors", hba, "Lipinski rule: <= 10"),
        ("Rotatable bonds", Descriptors.NumRotatableBonds(mol), "Flexibility"),
        ("Aromatic rings", rdMolDescriptors.CalcNumAromaticRings(mol), "Ring count"),
        ("Heavy atoms", mol.GetNumHeavyAtoms(), "Non-H atoms"),
        ("QED", round(QED.qed(mol), 3), "Drug-likeness, 0-1"),
        ("Lipinski violations", viol, "Ro5: 0-1 typical for drugs"),
        ("Murcko scaffold", MurckoScaffold.MurckoScaffoldSmiles(mol=mol), "Core ring system"),
    ], caption="Molecular descriptors")
    img = Draw.MolToImage(mol, size=(360, 260))
    st.image(img)


def full_smiles_vector(s):
    mol = Chem.MolFromSmiles(str(s).strip())
    if mol is None:
        st.error("Could not parse that SMILES.")
        return
    calc = get_mordred_calculator()
    df = calc.pandas([mol]).apply(pd.to_numeric, errors="coerce")
    df.insert(0, "SMILES", Chem.MolToSmiles(mol))
    offer_csv_download(df, "part2_molecule_features.csv")


def batch_smiles(df, col):
    keep, mols, failed = [], [], []
    for s in df[col].astype(str):
        s = s.strip()
        if not s or s.lower() == "nan":
            continue
        m = Chem.MolFromSmiles(s)
        if m is None:
            failed.append(s)
        else:
            mols.append(m)
            keep.append(Chem.MolToSmiles(m))
    if not mols:
        st.error("No valid SMILES found in that column.")
        return
    calc = get_mordred_calculator()
    feat = calc.pandas(mols).apply(pd.to_numeric, errors="coerce")
    feat.insert(0, "input_SMILES", keep)
    offer_excel_download(feat, "part2_molecule_FULL_features.xlsx", failed)


def chart_smiles(s):
    mol = Chem.MolFromSmiles(str(s).strip())
    if mol is None:
        st.error("Parse a valid SMILES first.")
        return
    mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
    hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    labels = ["MW/500", "LogP/5", "HBD/5", "HBA/10", "TPSA/140"]
    ratios = [max(0, mw / 500), max(0, logp / 5), hbd / 5, hba / 10, tpsa / 140]
    ang = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
    ang += ang[:1]
    vals = ratios + ratios[:1]
    fig = plt.figure(figsize=(6.6, 7.2))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.plot(ang, vals, color="#2563eb", linewidth=2)
    ax.fill(ang, vals, alpha=0.25, color="#2563eb")
    ax.plot(ang, [1] * len(ang), "--", color="#dc2626", linewidth=1.5,
            label="Rule-of-5 limit (inside the ring = passes)")
    ax.set_xticks(ang[:-1])
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_rlabel_position(36)
    ax.set_ylim(0, max(1.25, max(vals) * 1.1))
    fig.suptitle("Drug-likeness vs. rule-of-5", y=0.97, fontsize=14, fontweight="bold")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.07), frameon=False, fontsize=10)
    fig.subplots_adjust(top=0.84, bottom=0.15, left=0.10, right=0.90)
    st.pyplot(fig)
    st.download_button("Download 600-dpi PNG", data=save_fig_bytes(fig),
                        file_name="chart_molecule_druglikeness.png", mime="image/png")


# ---------------- UI ----------------
st.title("⚛️ Part 2 · Molecule Features")
st.caption("small-molecule SMILES  ·  RDKit + mordredcommunity  ·  1613 features")

st.markdown("**Single molecule — examples:**")
ex_cols = st.columns(3)
examples = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
}
if "p2_smiles" not in st.session_state:
    st.session_state.p2_smiles = "CC(=O)Oc1ccccc1C(=O)O"
for col, (label, value) in zip(ex_cols, examples.items()):
    if col.button(label):
        st.session_state.p2_smiles = value

smiles_input = st.text_area("Paste a SMILES string...", key="p2_smiles", height=70)

b1, b2, b3 = st.columns(3)
do_analyze = b1.button("Analyze", type="primary")
do_vector = b2.button("Full vector → CSV")
do_chart = b3.button("Charts")

if do_analyze:
    analyze_smiles(smiles_input)
if do_vector:
    full_smiles_vector(smiles_input)
if do_chart:
    chart_smiles(smiles_input)

st.markdown("---")
st.markdown("**Batch a dataset**")
st.caption("Upload a CSV / Excel / JSON with a column of SMILES (e.g. `smiles`). "
           "Returns one Excel row per molecule x 1613 descriptors.")

colname = st.text_input("Column name (blank = auto-detect)", key="p2_col")
uploaded = st.file_uploader("Upload list → full-feature Excel", type=["csv", "xlsx", "xls", "json"], key="p2_upload")

if uploaded is not None and st.button("Process batch", key="p2_run_batch"):
    df = read_uploaded_table(uploaded)
    if df is not None:
        col = colname.strip() or pick_column(df, SMILES_COLS, "SMILES")
        st.write(f"Featurizing column '{col}' with Mordred (1613 descriptors) ...")
        batch_smiles(df, col)
