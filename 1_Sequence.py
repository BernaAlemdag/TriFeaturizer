import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import peptides

from utils.shared import (
    render_feature_table, offer_csv_download, offer_excel_download,
    pick_column, read_uploaded_table, save_fig_bytes,
)

st.set_page_config(page_title="Sequence Features", page_icon="🧬", layout="wide")

plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 600, "savefig.bbox": "tight", "font.size": 11})

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
SEQ_COLS = ["sequence", "seq", "peptide", "protein", "fasta", "sequences"]

AA_TABLE = {
    'A': ('Alanine', 'Ala', 89.09, 6.00, 1.8, 'nonpolar'), 'R': ('Arginine', 'Arg', 174.20, 10.76, -4.5, 'positive'),
    'N': ('Asparagine', 'Asn', 132.12, 5.41, -3.5, 'polar'), 'D': ('Aspartate', 'Asp', 133.10, 2.77, -3.5, 'negative'),
    'C': ('Cysteine', 'Cys', 121.16, 5.07, 2.5, 'polar'), 'Q': ('Glutamine', 'Gln', 146.15, 5.65, -3.5, 'polar'),
    'E': ('Glutamate', 'Glu', 147.13, 3.22, -3.5, 'negative'), 'G': ('Glycine', 'Gly', 75.07, 5.97, -0.4, 'nonpolar'),
    'H': ('Histidine', 'His', 155.16, 7.59, -3.2, 'positive'), 'I': ('Isoleucine', 'Ile', 131.17, 6.02, 4.5, 'nonpolar'),
    'L': ('Leucine', 'Leu', 131.17, 5.98, 3.8, 'nonpolar'), 'K': ('Lysine', 'Lys', 146.19, 9.74, -3.9, 'positive'),
    'M': ('Methionine', 'Met', 149.21, 5.74, 1.9, 'nonpolar'), 'F': ('Phenylalanine', 'Phe', 165.19, 5.48, 2.8, 'nonpolar'),
    'P': ('Proline', 'Pro', 115.13, 6.30, -1.6, 'nonpolar'), 'S': ('Serine', 'Ser', 105.09, 5.68, -0.8, 'polar'),
    'T': ('Threonine', 'Thr', 119.12, 5.60, -0.7, 'polar'), 'W': ('Tryptophan', 'Trp', 204.23, 5.89, -0.9, 'nonpolar'),
    'Y': ('Tyrosine', 'Tyr', 181.19, 5.66, -1.3, 'polar'), 'V': ('Valine', 'Val', 117.15, 5.96, 4.2, 'nonpolar'),
}


def clean_sequence(text):
    body = "".join(l for l in str(text).strip().splitlines() if not l.strip().startswith(">"))
    return "".join(c for c in body.upper() if c.isalpha())


def featurize_sequence_full(text):
    clean = "".join(c for c in clean_sequence(text) if c in STANDARD_AA)
    if not clean:
        return None
    d = peptides.Peptide(clean).descriptors()
    pa = ProteinAnalysis(clean)
    d.update({"Length": len(clean), "MW": pa.molecular_weight(),
              "pI": pa.isoelectric_point(), "GRAVY": pa.gravy(),
              "Instability": pa.instability_index(), "Aromaticity": pa.aromaticity()})
    return d


def analyze_sequence(text):
    seq = clean_sequence(text)
    if not seq:
        st.warning("Please paste a sequence (raw or FASTA).")
        return
    if len(seq) == 1 and seq in AA_TABLE:
        name, tlc, mw, pI, kd, cls = AA_TABLE[seq]
        render_feature_table([
            ("Amino acid", f"{name} ({tlc}, {seq})", "Residue identity"),
            ("Molecular weight (Da)", mw, "Free amino-acid mass"),
            ("Isoelectric point (pI)", pI, "pH of zero net charge"),
            ("Hydropathy (Kyte-Doolittle)", kd, "+ hydrophobic / - hydrophilic"),
            ("Side-chain class", cls, "Polarity / charge"),
        ], caption=f"Single residue · {name}")
        return
    nonstd = sorted(set(seq) - STANDARD_AA)
    if nonstd:
        st.warning(f"Ignoring non-standard characters: {nonstd}")
    clean = "".join(c for c in seq if c in STANDARD_AA)
    if not clean:
        st.error("No standard amino acids found.")
        return
    pa = ProteinAnalysis(clean)
    ext = pa.molar_extinction_coefficient()
    rows = [
        ("Length (residues)", len(clean), "Chain length"),
        ("Molecular weight (Da)", round(pa.molecular_weight(), 2), "Peptide / protein mass"),
        ("Isoelectric point (pI)", round(pa.isoelectric_point(), 2), "pH of zero net charge"),
        ("Net charge @ pH 7.0", round(pa.charge_at_pH(7.0), 2), "Charge at neutral pH"),
        ("GRAVY", round(pa.gravy(), 3), "Avg hydropathy; + = hydrophobic"),
        ("Aromaticity", round(pa.aromaticity(), 3), "Fraction Phe+Trp+Tyr"),
        ("Instability index", round(pa.instability_index(), 2), "> 40 => likely unstable in vitro"),
        ("Extinction coeff (reduced)", ext[0], "M-1cm-1 @280 nm, Cys reduced"),
    ]
    pep = peptides.Peptide(clean)
    try:
        hm = round(pep.hydrophobic_moment(window=min(11, len(clean))), 3)
    except Exception:
        hm = "n/a"
    rows += [
        ("Aliphatic index", round(pep.aliphatic_index(), 2), "Thermostability proxy (A,V,I,L)"),
        ("Boman index", round(pep.boman(), 2), "Protein-binding potential (kcal/mol)"),
        ("Hydrophobic moment", hm, "Amphipathicity"),
    ]
    render_feature_table(rows, caption=f"Sequence features · {len(clean)} residues")


def full_sequence_vector(text):
    d = featurize_sequence_full(text)
    if d is None:
        st.error("Please paste a valid sequence first.")
        return
    offer_csv_download(pd.DataFrame([d]), "part1_sequence_features.csv")


def batch_sequences(df, col):
    rows, failed = [], []
    for v in df[col].astype(str):
        v = v.strip()
        if not v or v.lower() == "nan":
            continue
        d = featurize_sequence_full(v)
        if d is None:
            failed.append(v)
            continue
        rows.append({"input_sequence": v, **d})
    if not rows:
        st.error("No valid sequences found in that column.")
        return
    offer_excel_download(pd.DataFrame(rows), "part1_sequence_FULL_features.xlsx", failed)


def chart_sequence(text):
    clean = "".join(c for c in clean_sequence(text) if c in STANDARD_AA)
    if len(clean) < 2:
        st.warning("Need a peptide/protein (>=2 residues) to draw charts.")
        return
    order = "ACDEFGHIKLMNPQRSTVWY"
    counts = Counter(clean)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    ax[0].bar(list(order), [counts.get(a, 0) for a in order], color="#2563eb")
    ax[0].set_title("Amino-acid composition")
    ax[0].set_ylabel("Count")
    kd = {a: AA_TABLE[a][4] for a in order}
    w = min(9, len(clean))
    prof = [sum(kd[clean[i + j]] for j in range(w)) / w for i in range(len(clean) - w + 1)]
    ax[1].plot(range(1, len(prof) + 1), prof, color="#dc2626")
    ax[1].axhline(0, color="gray", lw=0.8)
    ax[1].set_title(f"Kyte-Doolittle hydropathy (window {w})")
    ax[1].set_xlabel("Residue position")
    ax[1].set_ylabel("Hydropathy")
    fig.suptitle(f"Sequence charts · {len(clean)} residues", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    st.pyplot(fig)
    st.download_button("Download 600-dpi PNG", data=save_fig_bytes(fig),
                        file_name="chart_sequence.png", mime="image/png")


# ---------------- UI ----------------
st.title("🧬 Part 1 · Sequence Features")
st.caption("residue / peptide / protein  ·  Biopython + peptides.py  ·  ≈108 features")

st.markdown("**Single sequence — examples:**")
ex_cols = st.columns(3)
examples = {
    "β-amyloid frag": "DAEFRHDSGYEVHHQK",
    "Insulin A-chain": "GIVEQCCTSICSLYQLENYCN",
    "Single residue (W)": "W",
}
if "p1_seq" not in st.session_state:
    st.session_state.p1_seq = "DAEFRHDSGYEVHHQK"
for col, (label, value) in zip(ex_cols, examples.items()):
    if col.button(label):
        st.session_state.p1_seq = value

seq_input = st.text_area(
    "Paste a residue, peptide, or protein sequence (raw or FASTA)...",
    key="p1_seq", height=100,
)

b1, b2, b3 = st.columns(3)
do_analyze = b1.button("Analyze", type="primary")
do_vector = b2.button("Full vector → CSV")
do_chart = b3.button("Charts")

if do_analyze:
    analyze_sequence(seq_input)
if do_vector:
    full_sequence_vector(seq_input)
if do_chart:
    chart_sequence(seq_input)

st.markdown("---")
st.markdown("**Batch a dataset**")
st.caption("Upload a CSV / Excel / JSON with a column of sequences (e.g. `sequence`). "
           "Returns one Excel row per sequence.")

colname = st.text_input("Column name (blank = auto-detect)", key="p1_col")
uploaded = st.file_uploader("Upload list → full-feature Excel", type=["csv", "xlsx", "xls", "json"], key="p1_upload")

if uploaded is not None and st.button("Process batch", key="p1_run_batch"):
    df = read_uploaded_table(uploaded)
    if df is not None:
        col = colname.strip() or pick_column(df, SEQ_COLS, "sequence")
        st.write(f"Featurizing column '{col}' ...")
        batch_sequences(df, col)
