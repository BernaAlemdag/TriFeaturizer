"""
Shared utilities for TriFeaturizer Streamlit app.
Equivalent to the "Shared utilities" cell in the original Colab notebook.
"""
import io
import pandas as pd
import streamlit as st


def render_feature_table(rows, caption=""):
    """rows = list of (feature, value, explanation) tuples."""
    df = pd.DataFrame(rows, columns=["Feature", "Value", "What it means"])
    if caption:
        st.markdown(f"#### {caption}")
    st.dataframe(df, hide_index=True, use_container_width=True)


def df_to_csv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8")


def df_to_excel_bytes(df: pd.DataFrame) -> bytes:
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        df.to_excel(writer, index=False)
    return buf.getvalue()


def offer_csv_download(df: pd.DataFrame, filename: str, label: str = "Download CSV"):
    st.success(f"Wrote {df.shape[0]} row(s) x {df.shape[1]} column(s)")
    st.download_button(
        label=label,
        data=df_to_csv_bytes(df),
        file_name=filename,
        mime="text/csv",
    )


def offer_excel_download(df: pd.DataFrame, filename: str, failed=None, label: str = "Download Excel"):
    st.success(f"Feature matrix: {df.shape[0]} row(s) x {df.shape[1]} column(s)")
    if failed:
        preview = ", ".join(map(str, failed[:3]))
        more = " ..." if len(failed) > 3 else ""
        st.warning(f"Skipped {len(failed)} invalid entr(y/ies): {preview}{more}")
    st.download_button(
        label=label,
        data=df_to_excel_bytes(df),
        file_name=filename,
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )


def pick_column(df: pd.DataFrame, candidates, label="input"):
    """Find the input column by common names; fall back to the first column."""
    low = {str(c).strip().lower(): c for c in df.columns}
    for n in candidates:
        if n in low:
            return low[n]
    st.info(f"No obvious {label} column found - using the first column '{df.columns[0]}'.")
    return df.columns[0]


def read_uploaded_table(uploaded_file):
    """Read an uploaded CSV / Excel / JSON file into a DataFrame."""
    if uploaded_file is None:
        return None
    name = uploaded_file.name
    low = name.lower()
    try:
        if low.endswith((".xlsx", ".xls")):
            df = pd.read_excel(uploaded_file)
        elif low.endswith(".json"):
            df = pd.read_json(uploaded_file)
        else:
            df = pd.read_csv(uploaded_file)
    except Exception as e:
        st.error(f"Could not read that file: {e}")
        return None
    st.info(f"Loaded {name}: {df.shape[0]} rows x {df.shape[1]} columns -> {list(df.columns)}")
    return df


def save_fig_bytes(fig, dpi=600) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()
