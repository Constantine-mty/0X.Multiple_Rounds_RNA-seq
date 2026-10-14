"""GSEA-specific GO BP classification and comparison summaries.

Biological category rules live in go_bp_classification.py. This module adds
GSEA-specific fields (absNES and Direction) and comparison aggregation.
"""

import pandas as pd

if __package__:
    from . import go_bp_classification as go_bp
else:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from gobp_pipeline import go_bp_classification as go_bp


# Compatibility names for existing callers of categorize_fine_en.
RULES = go_bp.RULES
THEME_RULES = go_bp.THEME_RULES
_CATEGORY_EN = go_bp._CATEGORY_EN
_BACKGROUND_CATEGORIES = go_bp._BACKGROUND_CATEGORIES
_HIERARCHY = go_bp._HIERARCHY
_HIERARCHY_AVAILABLE = go_bp._HIERARCHY_AVAILABLE
categorize_cn = go_bp.categorize_cn
categorize = go_bp.categorize
theme_of = go_bp.theme_of
is_signaling = go_bp.is_signaling
merged_category = go_bp.merged_category
_category_en = go_bp._category_en

COMPARISONS = []


def _try_load_hierarchy():
    """Refresh the shared GO hierarchy and the compatibility status fields."""
    go_bp._try_load_hierarchy()
    global _HIERARCHY, _HIERARCHY_AVAILABLE
    _HIERARCHY = go_bp._HIERARCHY
    _HIERARCHY_AVAILABLE = go_bp._HIERARCHY_AVAILABLE


def validate_gsea_result_columns(df: pd.DataFrame):
    """Check fields needed for GSEA-specific classification."""
    required = {"Description", "NES", "Comparison"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {', '.join(missing)}")
    if "Direction" not in df.columns and "GSEA_Direction" not in df.columns:
        raise ValueError("Input CSV needs either Direction or GSEA_Direction.")


def classify_gsea_results(df: pd.DataFrame) -> pd.DataFrame:
    """Apply GO BP classification, then add GSEA direction and absolute NES."""
    validate_gsea_result_columns(df)
    has_id = "ID" in df.columns
    classified = go_bp.classify_go_bp_table(df if has_id else df.assign(ID=None))
    if not has_id:
        classified = classified.drop(columns="ID")

    global COMPARISONS
    COMPARISONS = sorted(classified["Comparison"].dropna().unique().tolist())
    classified["absNES"] = classified["NES"].abs()
    if "Direction" not in classified.columns:
        classified["Direction"] = classified["GSEA_Direction"].map(
            {"Positive": "上调", "Negative": "下调"}
        )
    return classified


def load_gsea_results(path):
    """Read a GSEA CSV and classify it using the shared GO BP rules."""
    return classify_gsea_results(pd.read_csv(path))


def normalize_gsea_directions(df: pd.DataFrame) -> pd.Series:
    """Return analysis labels for Positive/Negative or Chinese up/down values."""
    column = "Direction" if "Direction" in df.columns else "GSEA_Direction"
    if column not in df.columns:
        raise ValueError("Input needs Direction or GSEA_Direction.")
    labels = {
        "positive": "上调", "up": "上调", "upregulated": "上调", "上调": "上调",
        "negative": "下调", "down": "下调", "downregulated": "下调", "下调": "下调",
    }
    normalized = df[column].astype("string").str.strip().str.casefold().map(labels)
    if normalized.isna().any():
        invalid = df.loc[normalized.isna(), column].drop_duplicates().tolist()
        raise ValueError(f"Unsupported {column} values: {invalid}")
    return normalized


def _empty_aggregate(category_col):
    return pd.DataFrame(columns=[
        category_col,
        "n",
        "median_absNES",
        "max_absNES",
        "min_padj",
        "top_term",
        "top_NES",
        "prominence",
    ])


def summarize_gsea_category(sub, category_col="Category"):
    """对单个 比较×方向 子集按指定类别列聚合。"""
    if category_col not in sub.columns:
        raise KeyError(f"Column not found: {category_col}")

    rows = []
    for cat, g in sub.groupby(category_col):
        rows.append({
            category_col: cat,
            "n": len(g),
            "median_absNES": round(g["absNES"].median(), 2),
            "max_absNES": round(g["absNES"].max(), 2),
            "min_padj": g["p.adjust"].min(),
            "top_term": g.sort_values("absNES", ascending=False).iloc[0]["Description"],
            "top_NES": round(g.sort_values("absNES", ascending=False).iloc[0]["NES"], 2),
            "prominence": round(len(g) * g["absNES"].median(), 2),
        })
    if not rows:
        return _empty_aggregate(category_col)
    agg = pd.DataFrame(rows).sort_values("prominence", ascending=False).reset_index(drop=True)
    return agg


def summarize_gsea_comparisons(df, category_col="Category"):
    """Summarize comparisons and directions within one Group only."""
    if category_col not in df.columns:
        raise KeyError(f"Column not found: {category_col}")
    if "Group" in df.columns and df["Group"].nunique(dropna=False) > 1:
        raise ValueError("Multiple Groups found; use summarize_gsea_groups().")

    directions = normalize_gsea_directions(df)
    source = df if "absNES" in df.columns else df.assign(absNES=df["NES"].abs())
    out = {}
    comps = sorted(df["Comparison"].dropna().unique().tolist())
    for comp in comps:
        for direction in ["上调", "下调"]:
            sub = source[(source["Comparison"] == comp) & (directions == direction)]
            out[(comp, direction)] = summarize_gsea_category(sub, category_col=category_col)
    return out


def summarize_gsea_groups(df, category_col="Category"):
    """Return comparison/direction statistics separately for each Group."""
    if "Group" not in df.columns:
        raise ValueError("Input is missing Group.")
    names = df["Group"].astype("string").str.strip()
    if names.isna().any() or names.eq("").any():
        raise ValueError("Group contains missing or blank values.")
    return {
        str(group): summarize_gsea_comparisons(group_df, category_col=category_col)
        for group, group_df in df.assign(Group=names).groupby("Group", sort=False)
    }


def track_gsea_themes(df, category_col="Category"):
    """5 主题 × 比较 的跨轮次追踪表，并报告 top term 对应的指定类别列。"""
    if category_col not in df.columns:
        raise KeyError(f"Column not found: {category_col}")

    records = []
    tdf = df[df["Theme"].notna()].copy()
    tdf["_Direction_zh"] = normalize_gsea_directions(tdf)
    comps = sorted(df["Comparison"].dropna().unique().tolist())
    for theme in [t for t, _ in THEME_RULES]:
        for comp in comps:
            sub = tdf[(tdf.Theme == theme) & (tdf.Comparison == comp)]
            if len(sub) == 0:
                records.append({"Theme": theme, "Comparison": comp, "Direction": "—",
                                "n": 0, "max_absNES": None, "top_term": "—",
                                "top_NES": None, "top_category": "—"})
            else:
                up = sub[sub["_Direction_zh"] == "上调"]
                dn = sub[sub["_Direction_zh"] == "下调"]
                if len(up) == 0:
                    strong, dlabel = dn, "下调"
                elif len(dn) == 0:
                    strong, dlabel = up, "上调"
                else:
                    if up["absNES"].max() >= dn["absNES"].max():
                        strong, dlabel = up, "上调"
                    else:
                        strong, dlabel = dn, "下调"
                top = strong.sort_values("absNES", ascending=False).iloc[0]
                records.append({"Theme": theme, "Comparison": comp, "Direction": dlabel,
                                "n": len(sub), "max_absNES": round(sub["absNES"].max(), 2),
                                "top_term": top["Description"],
                                "top_NES": round(top["NES"], 2),
                                "top_category": top[category_col]})
    return pd.DataFrame(records)


def main(argv=None):
    """Print a compact classification check for one GSEA CSV."""
    import sys
    args = sys.argv[1:] if argv is None else argv
    path = args[0] if args else "02.GSEA_H_fdr0.05.csv"
    df = load_gsea_results(path)
    print(f"Hierarchy available: {_HIERARCHY_AVAILABLE}")
    print(f"Total rows: {len(df)}, Unique descriptions: {df['Description'].nunique()}")
    vc = df["Category"].value_counts()
    print(f"\nCategory distribution (n={len(vc)} cats):")
    print(vc.to_string())
    other_pct = vc.get("Other/Uncategorized", 0) / len(df) * 100
    bg_total = df["Background"].sum()
    bg_pct = bg_total / len(df) * 100
    print(f"\nOther/Uncategorized: {other_pct:.1f}%   背景集总计: {bg_pct:.1f}% ({bg_total} rows)")
    print("\nOther/Uncategorized sample (first 30):")
    others = df[df.Category == "Other/Uncategorized"].drop_duplicates("Description")
    for _, r in others.head(30).iterrows():
        print(f"  {r['ID']}  {r['Description']}")


# Backward-compatible API names for existing notebooks and scripts.
validate_gsea_columns = validate_gsea_result_columns
classify_gsea = classify_gsea_results
load = load_gsea_results
aggregate_group = summarize_gsea_category
all_aggregates = summarize_gsea_comparisons
theme_tracking = track_gsea_themes


if __name__ == "__main__":
    main()
