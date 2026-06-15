# 🧪 TriFeaturizer

**A multi-domain molecular feature extractor.** Turn protein/peptide sequences, small-molecule SMILES, and inorganic materials (chemical formulas or CIF crystal structures) into complete, machine-learning-ready feature sets one item at a time, or for a whole uploaded dataset at once.

![Python](https://img.shields.io/badge/python-3.10%2B-blue)

### ▶ [Run it live in Google Colab](https://colab.research.google.com/github/USERNAME/REPO/blob/main/TriFeaturizer.v2.ipynb)

> The interactive interface (buttons, file upload, charts) only works in a **live kernel** — open the notebook in **Colab** to use it. GitHub shows the *code* but cannot run widgets, so the buttons won't be clickable there. For a populated, non-interactive view on GitHub, see [`TriFeaturizer_preview.ipynb`](TriFeaturizer_preview.ipynb), which has the example outputs (tables, molecule image, 600-dpi charts) saved in.

> Replace `USERNAME/REPO` in the Colab links and badge with your own GitHub path after you upload.

---

## Overview

TriFeaturizer is a single Google Colab notebook that unifies three different chemical "languages" behind one consistent interface. Whatever you feed it, you get back two things: a short **explained table** of the key descriptors (so you understand the numbers), and the **full descriptor vector** ready for modelling. For datasets, it produces a single Excel workbook with one row per item and one column per feature.

Within each domain, two engines work together one provides the readable values, the other emits the large vector:

| Part | Input | Engines | Features produced |
|------|-------|---------|-------------------|
| **1 · Sequence** | residue / peptide / protein | Biopython + peptides.py | ≈ 108 |
| **2 · Molecule** | SMILES | RDKit + mordredcommunity | 1613 |
| **3 · Material** | chemical formula **and/or** CIF | pymatgen + matminer + DScribe | 142 (+ SOAP for CIF) |

---

## Features

- **Three input domains, one notebook** Sequences, molecules, and inorganic materials.
- **Explained + full views**: A curated table for humans, and the complete vector for models.
- **Batch mode → Excel** : Upload a CSV / Excel / JSON list and download a full feature matrix as `.xlsx`.
- **Charts** : Amino-acid composition & Kyte–Doolittle hydropathy, drug-likeness vs. rule-of-five radar, element fractions & lattice parameters all exported as **600-dpi PNGs**.
- **CIF crystal support** : Parse a structure, get density/symmetry features and an optional SOAP fingerprint.
- **Per-section installs** : A sequence-only user never pulls the heavy materials stack.
- **Companion reference** : A Word document listing every feature name and meaning.

---

## Quickstart

### Run in Colab (recommended)

1. Click the **Open in Colab** badge above (after setting your repo path).
2. Run the **Shared utilities** cell once.
3. For each part you need, run its `pip install` cell, then its code cell.
4. Use the on-screen buttons paste an item, or upload a dataset.

### Run locally

```bash
pip install biopython peptides rdkit mordredcommunity \
            pymatgen matminer dscribe openpyxl ipywidgets pandas matplotlib
jupyter notebook TriFeaturizer.v2.ipynb
```

> Note: the **Upload** and **Download** buttons use Colab APIs. Locally, load your data into a pandas DataFrame and call the `batch_*` / `featurize_*` functions directly.

---

## Usage

### A · Explore one item

Paste a single sequence, SMILES, or formula, then choose:

- **Analyze** An explained table of the key descriptors.
- **Charts** Quick plots (saved at 600 dpi).
- **Full vector → CSV** The complete descriptor set for that one item.

### B · Process a dataset into Excel

Upload a **CSV, Excel (.xlsx), or JSON** file containing a column of inputs. The relevant column is auto-detected by name (`smiles`, `sequence`, `formula`, …); you can also type the column name. TriFeaturizer featurizes every row and returns one workbook.

**Input**  One row per item; other columns (IDs, labels, targets) are ignored and can be re-joined later:

```csv
id,smiles
mol_1,CC(=O)Oc1ccccc1C(=O)O
mol_2,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
```

**Output** First column echoes the input, every other column is a descriptor:

| Workbook | Rows × Columns |
|----------|----------------|
| `part1_sequence_FULL_features.xlsx` | N sequences × (1 + ~108) |
| `part2_molecule_FULL_features.xlsx` | N molecules × (1 + 1613) |
| `part3_composition_FULL_features.xlsx` | N formulas × (1 + 142) |

Rows that can't be parsed are reported and skipped (they don't stop the run); descriptors that are undefined for an otherwise-valid item appear as blanks. The file loads straight into pandas / scikit-learn.

---

## Feature counts

| Domain | Full-vector size | Breakdown |
|--------|------------------|-----------|
| Sequence | ≈ 108 | 102 peptides.py QSAR descriptors + 6 Biopython properties |
| Molecule | 1613 | Mordred 2-D descriptors |
| Composition | 142 | 132 Magpie (22 element properties × 6 statistics) + 4 valence-orbital fractions + 6 stoichiometry norms |
| CIF structure | 8 + SOAP | 3 density + 5 symmetry features, plus an optional SOAP block (≈ 390–1500, growing with the number of distinct elements) |

Every individual feature name and meaning is documented in **[`docs/TriFeaturizer_Feature_Reference.docx`](docs/TriFeaturizer_Feature_Reference.docx)**.

---

## Repository structure

```
TriFeaturizer/
├── TriFeaturizer_v2 (1).ipynb            # the interactive notebook (run in Colab)
├── TriFeaturizer_preview.ipynb         # static preview with example outputs (renders on GitHub)
├── README.md
├── docs/
│   ├── TriFeaturizer_Manual.pdf   # full feature dictionary
└── examples/
    ├── sequences.csv
    ├── molecules.csv
    └── formulas.csv
```
---

## Built with

TriFeaturizer is a thin, user-friendly layer over established open-source libraries. Please cite the underlying tools when you publish results.

| Library | Used for |
|---------|----------|
| [Biopython](https://biopython.org/) | sequence physico-chemical properties |
| [peptides.py](https://github.com/althonos/peptides.py) | QSAR sequence descriptors |
| [RDKit](https://www.rdkit.org/) | molecule descriptors & depictions |
| [Mordred (community fork)](https://github.com/JacksonBurns/mordred-community) | full 2-D molecular descriptor vector |
| [pymatgen](https://pymatgen.org/) | composition & crystal-structure parsing |
| [matminer](https://hackingmaterials.lbl.gov/matminer/) | Magpie / valence / density / symmetry featurizers |
| [DScribe](https://singroup.github.io/dscribe/) | SOAP structural fingerprint |

---

## Troubleshooting

**`pip's dependency resolver … google-colab 1.0.0 requires requests==2.32.4, but you have requests 2.34.2`**

This is a **harmless warning, not an error**. A materials-stack dependency upgrades `requests`, which no longer matches the version `google-colab` declares. Feature extraction is unaffected. To silence it, run `!pip -q install "requests==2.32.4"` and then **Runtime ▸ Restart session**, or install with the pin: `!pip -q install matminer dscribe "requests==2.32.4"`.

**matminer / pymatgen install takes a couple of minutes** that's expected; it pulls a large dependency chain (pymatgen + ASE).

**Materials warning**  matminer is built for *inorganic* compositions and crystals. Rule of thumb: organic / connectivity matters → Part 2 (SMILES); oxide, alloy, crystal → Part 3.

## Citation

If TriFeaturizer helps your work, please cite this repository and the underlying libraries listed above.
