# TriFeaturizer 🧪

Multi-domain molecular feature extractor — web app version, converted from the
original Colab notebook (`TriFeaturizer_v2.ipynb`).

Turn a **protein/peptide sequence**, a **small-molecule SMILES**, or an
**inorganic material** (chemical formula or CIF crystal structure) into a
machine-learning-ready feature set — one item at a time, or as a batch
Excel export.

## Repository structure

```
trifeaturizer/
├── App.py                  # Entry point / landing page
├── requirements.txt        # All Python dependencies
├── pages/                  # Each file here = one sidebar page (Streamlit convention)
│   ├── 1_Sequence.py        # Part 1: Biopython + peptides.py  (~108 features)
│   ├── 2_Molecule.py        # Part 2: RDKit + mordredcommunity (1613 features)
│   └── 3_Material.py        # Part 3: pymatgen + matminer + DScribe (142 + SOAP)
└── utils/
    └── shared.py            # Table rendering, CSV/Excel download helpers
```

## Run locally

```bash
pip install -r requirements.txt
streamlit run App.py
```

This opens a browser tab at `http://localhost:8501`. Editing any file and
saving auto-reloads the page.

## Deploy

**Option A — Render**
1. Push this folder to a GitHub repo.
2. On Render, create a new "Web Service", connect the repo.
3. Build command: `pip install -r requirements.txt`
4. Start command: `streamlit run App.py --server.port $PORT --server.address 0.0.0.0`

**Option B — Streamlit Community Cloud**
1. Push this folder to a GitHub repo.
2. Go to share.streamlit.io, connect the repo, point it at `App.py`.
3. Done — no start command needed, it's built for this.

## Note on install size

`pymatgen`, `matminer`, and `dscribe` (used in Part 3 / Material) are heavier,
partly-compiled packages. Expect the first build/deploy to take a few minutes
longer than a typical lightweight Streamlit app.
