import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

from pymatgen.core import Composition, Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from matminer.featurizers.composition import ElementProperty, ValenceOrbital, Stoichiometry

from utils.shared import (
    render_feature_table, offer_csv_download, offer_excel_download,
    pick_column, read_uploaded_table, save_fig_bytes,
)

st.set_page_config(page_title="Material Features", page_icon="🪨", layout="wide")

plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 600, "savefig.bbox": "tight", "font.size": 11})

FORMULA_COLS = ["formula", "composition", "compound", "formulae", "compounds", "material"]

CURATED_MAGPIE = [
    ("MagpieData mean Number", "Mean atomic number"),
    ("MagpieData mean AtomicWeight", "Mean atomic weight"),
    ("MagpieData mean Electronegativity", "Mean Pauling electronegativity"),
    ("MagpieData mean CovalentRadius", "Mean covalent radius (pm)"),
    ("MagpieData mean NValence", "Mean valence-electron count"),
    ("MagpieData avg_dev Electronegativity", "Spread of electronegativity"),
]


@st.cache_resource
def get_featurizers():
    ep = ElementProperty.from_preset("magpie")
    vo = ValenceOrbital(props=["frac"])
    st_ = Stoichiometry()
    return ep, vo, st_


def featurize_composition_full(formula):
    try:
        comp = Composition(formula)
    except Exception:
        return None
    ep, vo, st_ = get_featurizers()
    d = {}
    for f in (ep, vo, st_):
        d.update(dict(zip(f.feature_labels(), f.featurize(comp))))
    return d


# ---- 3a composition ----
def analyze_composition(formula):
    try:
        comp = Composition(formula)
    except Exception as e:
        st.error(f"Could not parse formula: {e}")
        return
    ep, vo, _ = get_featurizers()
    ep_d = dict(zip(ep.feature_labels(), ep.featurize(comp)))
    vo_d = dict(zip(vo.feature_labels(), vo.featurize(comp)))
    rows = [
        ("Reduced formula", comp.reduced_formula, "Normalised composition"),
        ("Elements", ", ".join(sorted(e.symbol for e in comp.elements)), "Constituents"),
        ("Atoms / formula unit", round(comp.num_atoms, 3), "Total atoms"),
        ("Avg atomic mass", round(float(comp.weight) / comp.num_atoms, 3), "amu per atom"),
    ]
    for lbl, meaning in CURATED_MAGPIE:
        if lbl in ep_d:
            rows.append((lbl.replace("MagpieData ", "Magpie "), round(ep_d[lbl], 3), meaning))
    for lbl in vo.feature_labels():
        rows.append((lbl, round(vo_d[lbl], 3), "Fraction of valence e- in this orbital"))
    render_feature_table(rows, caption=f"Composition features · {comp.reduced_formula}")


def full_composition_vector(formula):
    d = featurize_composition_full(formula)
    if d is None:
        st.error("Could not parse that formula.")
        return
    df = pd.DataFrame([d])
    df.insert(0, "formula", Composition(formula).reduced_formula)
    offer_csv_download(df, "part3_composition_features.csv")


def batch_compositions(df, col):
    rows, failed = [], []
    for v in df[col].astype(str):
        v = v.strip()
        if not v or v.lower() == "nan":
            continue
        d = featurize_composition_full(v)
        if d is None:
            failed.append(v)
            continue
        rows.append({"input_formula": v, **d})
    if not rows:
        st.error("No valid formulas found in that column.")
        return
    offer_excel_download(pd.DataFrame(rows), "part3_composition_FULL_features.xlsx", failed)


def chart_composition(formula):
    try:
        comp = Composition(formula)
    except Exception as e:
        st.error(f"Could not parse formula: {e}")
        return
    fr = comp.fractional_composition.get_el_amt_dict()
    els = list(fr.keys())
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(els, [fr[e] for e in els], color="#0891b2")
    ax.set_title(f"Atomic fraction · {comp.reduced_formula}")
    ax.set_ylabel("Fraction")
    ax.set_ylim(0, 1)
    fig.tight_layout()
    st.pyplot(fig)
    st.download_button("Download 600-dpi PNG", data=save_fig_bytes(fig),
                        file_name="chart_composition.png", mime="image/png")


# ---- 3b structure / CIF ----
@st.cache_resource
def get_nacl_structure():
    return Structure.from_spacegroup(
        "Fm-3m", Lattice.cubic(5.64), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]]
    )


def analyze_structure(struct):
    sga = SpacegroupAnalyzer(struct)
    a, b, c = struct.lattice.abc
    render_feature_table([
        ("Formula", struct.composition.reduced_formula, "Structure composition"),
        ("Sites in cell", len(struct), "Atoms in the unit cell"),
        ("Density (g/cm³)", round(float(struct.density), 3), "Mass / volume"),
        ("Cell volume (Å³)", round(struct.volume, 2), "Unit-cell volume"),
        ("Volume / atom (Å³)", round(struct.volume / len(struct), 2), "Packing measure"),
        ("Lattice a,b,c (Å)", f"{a:.3f}, {b:.3f}, {c:.3f}", "Cell edge lengths"),
        ("Space group", f"{sga.get_space_group_symbol()} (#{sga.get_space_group_number()})", "Symmetry"),
        ("Crystal system", sga.get_crystal_system(), "Lattice family"),
    ], caption="Crystal-structure features")


def full_structure_vector(struct, with_soap=False):
    from matminer.featurizers.structure import DensityFeatures, GlobalSymmetryFeatures
    data = {}
    for f in (DensityFeatures(), GlobalSymmetryFeatures()):
        try:
            data.update(dict(zip(f.feature_labels(), f.featurize(struct))))
        except Exception as e:
            st.warning(f"Skipped {f.__class__.__name__}: {e}")
    if with_soap:
        from dscribe.descriptors import SOAP
        from pymatgen.io.ase import AseAtomsAdaptor
        atoms = AseAtomsAdaptor.get_atoms(struct)
        species = sorted(set(atoms.get_chemical_symbols()))
        soap = SOAP(species=species, r_cut=5.0, n_max=6, l_max=4, periodic=True)
        vec = soap.create(atoms).mean(axis=0)
        data.update({f"SOAP_{i}": v for i, v in enumerate(vec)})
    df = pd.DataFrame([data])
    df.insert(0, "formula", struct.composition.reduced_formula)
    offer_csv_download(df, "part3_structure_features.csv")


def chart_structure(struct):
    a, b, c = struct.lattice.abc
    fr = struct.composition.fractional_composition.get_el_amt_dict()
    els = list(fr.keys())
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].bar(els, [fr[e] for e in els], color="#0891b2")
    ax[0].set_title(f"Atomic fraction · {struct.composition.reduced_formula}")
    ax[0].set_ylabel("Fraction")
    ax[0].set_ylim(0, 1)
    ax[1].bar(["a", "b", "c"], [a, b, c], color="#7c3aed")
    ax[1].set_title("Lattice parameters (Å)")
    fig.tight_layout()
    st.pyplot(fig)
    st.download_button("Download 600-dpi PNG", data=save_fig_bytes(fig),
                        file_name="chart_structure.png", mime="image/png")


def load_structure(cif_text):
    return Structure.from_str(cif_text, fmt="cif")


# ==================== UI ====================
st.title("🪨 Part 3 · Material Features")
st.caption("chemical formula and/or CIF crystal structure  ·  pymatgen + matminer + DScribe  ·  142 features (+ optional SOAP)")
st.info("matminer is built for **inorganic** compositions and crystals. "
        "Rule of thumb: organic / connectivity matters -> use Part 2 (SMILES); oxide, alloy, crystal -> use this page.")

tab_a, tab_b = st.tabs(["3a · Composition (formula only)", "3b · Crystal structure (CIF)"])

# ---- Tab 3a: composition ----
with tab_a:
    st.markdown("**Examples:**")
    ex_cols = st.columns(3)
    examples = {"LiFePO₄": "LiFePO4", "Fe₂O₃": "Fe2O3", "BaTiO₃": "BaTiO3"}
    if "p3c_formula" not in st.session_state:
        st.session_state.p3c_formula = "LiFePO4"
    for col, (label, value) in zip(ex_cols, examples.items()):
        if col.button(label, key=f"ex_{label}"):
            st.session_state.p3c_formula = value

    formula_input = st.text_input("Chemical formula, e.g. Fe2O3", key="p3c_formula")

    b1, b2, b3 = st.columns(3)
    if b1.button("Analyze", type="primary", key="p3c_a"):
        analyze_composition(formula_input)
    if b2.button("Full vector → CSV", key="p3c_f"):
        full_composition_vector(formula_input)
    if b3.button("Charts", key="p3c_c"):
        chart_composition(formula_input)

    st.markdown("---")
    st.markdown("**Batch a dataset**")
    st.caption("Upload a CSV / Excel / JSON with a column of formulas (e.g. `formula`). "
               "Returns one Excel row per compound x 142 descriptors.")
    colname = st.text_input("Column name (blank = auto-detect)", key="p3c_col")
    uploaded = st.file_uploader("Upload list → full-feature Excel", type=["csv", "xlsx", "xls", "json"], key="p3c_upload")
    if uploaded is not None and st.button("Process batch", key="p3c_run_batch"):
        df = read_uploaded_table(uploaded)
        if df is not None:
            col = colname.strip() or pick_column(df, FORMULA_COLS, "formula")
            st.write(f"Featurizing column '{col}' (142 composition descriptors) ...")
            batch_compositions(df, col)

# ---- Tab 3b: structure / CIF ----
with tab_b:
    if "p3s_struct" not in st.session_state:
        st.session_state.p3s_struct = None
    if "p3s_cif_text" not in st.session_state:
        st.session_state.p3s_cif_text = ""

    c1, c2 = st.columns(2)
    if c1.button("Load NaCl example"):
        from pymatgen.io.cif import CifWriter
        nacl = get_nacl_structure()
        st.session_state.p3s_struct = nacl
        st.session_state.p3s_cif_text = str(CifWriter(nacl))
    uploaded_cif = c2.file_uploader("Upload CIF", type=["cif"], key="p3s_upload_cif")
    if uploaded_cif is not None:
        st.session_state.p3s_cif_text = uploaded_cif.read().decode("utf-8", errors="ignore")

    cif_text = st.text_area(
        "Paste CIF text here, or click \u201cLoad NaCl example\u201d / upload a CIF...",
        key="p3s_cif_text", height=160,
    )
    soap_checkbox = st.checkbox("Include SOAP (DScribe) in CSV", value=False, key="p3s_soap")

    def _resolve_structure():
        if cif_text.strip():
            try:
                return load_structure(cif_text)
            except Exception as e:
                st.error(f"Could not parse CIF: {e}")
                return None
        return st.session_state.p3s_struct

    b1, b2, b3 = st.columns(3)
    if b1.button("Analyze structure", type="primary", key="p3s_a"):
        s = _resolve_structure()
        if s is None:
            st.warning("Paste a CIF, upload one, or load the example.")
        else:
            analyze_structure(s)
    if b2.button("Full vector → CSV", key="p3s_f"):
        s = _resolve_structure()
        if s is None:
            st.warning("Paste a CIF, upload one, or load the example.")
        else:
            full_structure_vector(s, with_soap=soap_checkbox)
    if b3.button("Charts", key="p3s_c"):
        s = _resolve_structure()
        if s is None:
            st.warning("Paste a CIF, upload one, or load the example.")
        else:
            chart_structure(s)
