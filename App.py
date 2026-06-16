import streamlit as st

st.set_page_config(
    page_title="TriFeaturizer",
    page_icon="🧪",
    layout="wide",
)

st.title("🧪 TriFeaturizer")
st.subheader("Multi-domain Molecular Feature Extractor")

st.markdown(
    """
Turn a **protein/peptide sequence**, a **small-molecule SMILES**, or an
**inorganic material** (chemical formula or CIF crystal structure) into a
complete, machine-learning-ready feature set — either one item at a time,
or for a whole uploaded dataset at once.

Use the sidebar to pick a tool:
"""
)

col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("### 1 · Sequence")
    st.markdown(
        "Biopython + peptides.py\n\n"
        "**~108 features**\n\n"
        "Residue / peptide / protein sequences"
    )

with col2:
    st.markdown("### 2 · Molecule")
    st.markdown(
        "RDKit + mordredcommunity\n\n"
        "**1613 features**\n\n"
        "Small-molecule SMILES strings"
    )

with col3:
    st.markdown("### 3 · Material")
    st.markdown(
        "pymatgen + matminer + DScribe\n\n"
        "**142 features** (+ optional SOAP)\n\n"
        "Chemical formula and/or CIF crystal structure"
    )

st.markdown("---")
st.markdown(
    """
**Two ways to use each tool**

- **Explore one item.** Paste a single sequence, SMILES, or formula and click
  **Analyze** for the explained table, or **Full vector → CSV** for the
  complete descriptor set of that one item.
- **Process a whole dataset → Excel.** Upload a CSV, Excel, or JSON file
  containing a column of sequences (or SMILES, or formulas). TriFeaturizer
  computes the complete feature vector for every row and returns a single
  `.xlsx` workbook — one row per item, one column per descriptor.
"""
)
