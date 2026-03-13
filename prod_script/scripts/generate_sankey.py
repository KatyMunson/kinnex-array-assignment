#!/usr/bin/env python3
"""
generate_sankey.py

Generates an interactive Sankey diagram of read flow through the pipeline:
S-reads → Skera → assigned/unassigned → libraries → FLNC.
Output is a self-contained HTML file using Plotly.

Usage: python generate_sankey.py <results_qc_dir> <output_html>

Author:  KM
Created: 2026-02
"""

import sys
import json
import pandas as pd
from pathlib import Path

qc_dir  = Path(sys.argv[1])
out_html = Path(sys.argv[2])

# ── Load QC tables ─────────────────────────────────────────────────────────────
df_skera   = pd.read_csv(qc_dir / "qc_skera.tsv",       sep="\t")
df_split   = pd.read_csv(qc_dir / "qc_split_skera.tsv", sep="\t")
df_refine  = pd.read_csv(qc_dir / "qc_refine.tsv",      sep="\t")

# ── Build node and link lists ──────────────────────────────────────────────────
node_labels = []
node_colors = []
links_source = []
links_target = []
links_value  = []
links_color  = []

node_index = {}

def get_or_add_node(label, color):
    if label not in node_index:
        node_index[label] = len(node_labels)
        node_labels.append(label)
        node_colors.append(color)
    return node_index[label]

# Color palette
SAMPLE_COLOR    = "#2E86AB"
UNASSIGNED_COLOR = "#E84855"
LIBRARY_COLORS  = [
    "#3BB273", "#F4A261", "#9B5DE5", "#F15BB5",
    "#00BBF9", "#FEE440", "#00F5D4", "#FB8500",
    "#8338EC", "#06D6A0", "#FFB703", "#EF476F",
    "#118AB2", "#073B4C", "#264653", "#2A9D8F",
]
FLNC_COLOR      = "#457B9D"
LOST_COLOR      = "#BBBBBB"

# Distinct colors for bc01-bc12
BC_COLORS = {
    "bc01": "#E63946", "bc02": "#F4A261", "bc03": "#2A9D8F",
    "bc04": "#E9C46A", "bc05": "#264653", "bc06": "#9B5DE5",
    "bc07": "#F15BB5", "bc08": "#00BBF9", "bc09": "#06D6A0",
    "bc10": "#FB8500", "bc11": "#8338EC", "bc12": "#FF006E",
}

# Per-sample library color assignment (consistent across samples)
all_libs = sorted(df_split[df_split["library"] != "unassigned"]["library"].unique())
lib_color_map = {lib: LIBRARY_COLORS[i % len(LIBRARY_COLORS)] for i, lib in enumerate(all_libs)}

samples = df_skera["sample"].tolist()

for s in samples:
    sk_row = df_skera[df_skera["sample"] == s]
    if sk_row.empty:
        continue
    sreads = int(sk_row.iloc[0]["segmented_reads"])

    sample_node = get_or_add_node(f"{s}\nS-reads", SAMPLE_COLOR)

    sp = df_split[df_split["sample"] == s]
    sp_unassigned  = sp[sp["library"] == "unassigned"]
    sp_highconf    = sp[(sp["library"] != "unassigned") & (sp["confidence"] == "highconf")]
    sp_lowconf     = sp[(sp["library"] != "unassigned") & (sp["confidence"] == "lowconf")]

    assigned_total   = int(sp_highconf["read_count"].sum()) + int(sp_lowconf["read_count"].sum())
    unassigned_total = int(sp_unassigned["read_count"].sum()) if len(sp_unassigned) else 0
    lima_dropout     = sreads - assigned_total - unassigned_total

    if unassigned_total > 0:
        unassigned_node = get_or_add_node(f"{s}\nUnassigned", UNASSIGNED_COLOR)
        links_source.append(sample_node)
        links_target.append(unassigned_node)
        links_value.append(unassigned_total)
        links_color.append("rgba(232,72,85,0.3)")

    if lima_dropout > 0:
        dropout_node = get_or_add_node(f"{s}\nLima1 dropout", LOST_COLOR)
        links_source.append(sample_node)
        links_target.append(dropout_node)
        links_value.append(lima_dropout)
        links_color.append("rgba(187,187,187,0.3)")

    # Per-library, per-confidence flow
    all_libs_in_sample = sorted(set(sp[sp["library"] != "unassigned"]["library"].unique()))
    for lib in all_libs_in_sample:
        lcolor   = lib_color_map.get(lib, "#999999")
        lcolor_t = lcolor + "66"

        for conf_label, conf_alpha in [
            ("highconf", "99"),
            ("lowconf",  "55"),
        ]:
            lib_reads_rows = sp[(sp["library"] == lib) & (sp["confidence"] == conf_label)]
            if lib_reads_rows.empty:
                continue
            lib_reads = int(lib_reads_rows["read_count"].sum())
            if lib_reads == 0:
                continue

            lib_node = get_or_add_node(f"{lib}\n{conf_label}", lcolor)
            links_source.append(sample_node)
            links_target.append(lib_node)
            links_value.append(lib_reads)
            links_color.append(lcolor + "66")

            # FLNC for this sample + library + confidence
            rf = df_refine[
                (df_refine["sample"] == s) &
                (df_refine["library"] == lib) &
                (df_refine["confidence"] == conf_label)
            ]
            flnc_total = int(rf["total_flnc"].sum()) if len(rf) else 0
            flnc_lost  = lib_reads - flnc_total

            if flnc_total > 0:
                flnc_node = get_or_add_node(f"{lib}\n{conf_label}\nFLNC", FLNC_COLOR)
                links_source.append(lib_node)
                links_target.append(flnc_node)
                links_value.append(flnc_total)
                links_color.append("rgba(69,123,157,0.4)")

                # Per-barcode breakdown
                if len(rf):
                    for bc in [f"bc{i:02d}" for i in range(1, 13)]:
                        if bc not in rf.columns:
                            continue
                        bc_count = int(rf[bc].sum())
                        if bc_count > 0:
                            bc_color   = BC_COLORS.get(bc, "#999999")
                            bc_color_t = bc_color + "55"
                            bc_node = get_or_add_node(f"{lib}\n{conf_label}\n{bc}", bc_color)
                            links_source.append(flnc_node)
                            links_target.append(bc_node)
                            links_value.append(bc_count)
                            links_color.append(bc_color_t)

            if flnc_lost > 0:
                lost_node = get_or_add_node(f"{lib}\n{conf_label}\nNon-FLNC", LOST_COLOR)
                links_source.append(lib_node)
                links_target.append(lost_node)
                links_value.append(flnc_lost)
                links_color.append("rgba(187,187,187,0.3)")

# ── Build Plotly figure data as JSON ──────────────────────────────────────────
sankey_data = {
    "type": "sankey",
    "orientation": "h",
    "node": {
        "pad": 20,
        "thickness": 24,
        "line": {"color": "rgba(0,0,0,0.2)", "width": 0.5},
        "label": node_labels,
        "color": node_colors,
        "hovertemplate": "%{label}<br>%{value:,} reads<extra></extra>",
    },
    "link": {
        "source": links_source,
        "target": links_target,
        "value":  links_value,
        "color":  links_color,
        "hovertemplate": "%{source.label} → %{target.label}<br>%{value:,} reads<extra></extra>",
    }
}

layout = {
    "title": {
        "text": "Iso-Seq Pipeline Read Flow",
        "font": {"size": 22, "family": "IBM Plex Mono, monospace", "color": "#1a1a2e"},
        "x": 0.5,
        "xanchor": "center",
    },
    "font": {"family": "IBM Plex Mono, monospace", "size": 11, "color": "#333"},
    "paper_bgcolor": "#f8f9fb",
    "plot_bgcolor":  "#f8f9fb",
    "margin": {"l": 20, "r": 20, "t": 80, "b": 20},
    "height": max(600, len(samples) * 80),
}

fig_json = json.dumps({"data": [sankey_data], "layout": layout})

# ── Write HTML ────────────────────────────────────────────────────────────────
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Iso-Seq Pipeline Read Flow</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600&display=swap" rel="stylesheet">
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    background: #f8f9fb;
    font-family: 'IBM Plex Sans', sans-serif;
    color: #1a1a2e;
    min-height: 100vh;
    padding: 2rem;
  }}
  header {{
    margin-bottom: 1.5rem;
    border-left: 4px solid #2E86AB;
    padding-left: 1rem;
  }}
  header h1 {{
    font-family: 'IBM Plex Mono', monospace;
    font-size: 1.4rem;
    font-weight: 600;
    letter-spacing: -0.02em;
    color: #1a1a2e;
  }}
  header p {{
    font-size: 0.85rem;
    color: #666;
    margin-top: 0.25rem;
  }}
  .legend {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.75rem;
    margin-bottom: 1.5rem;
    padding: 1rem 1.25rem;
    background: white;
    border-radius: 6px;
    border: 1px solid #e5e7eb;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05);
  }}
  .legend-item {{
    display: flex;
    align-items: center;
    gap: 0.4rem;
    font-size: 0.78rem;
    font-family: 'IBM Plex Mono', monospace;
    color: #444;
  }}
  .legend-dot {{
    width: 12px;
    height: 12px;
    border-radius: 2px;
    flex-shrink: 0;
  }}
  #sankey {{
    background: white;
    border-radius: 6px;
    border: 1px solid #e5e7eb;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    padding: 1rem;
  }}
</style>
</head>
<body>
<header>
  <h1>Iso-Seq Pipeline — Read Flow</h1>
  <p>S-reads &rarr; ZMW assignment &rarr; libraries &rarr; FLNC</p>
</header>
<div class="legend">
  <div class="legend-item"><div class="legend-dot" style="background:#2E86AB"></div>Sample S-reads</div>
  <div class="legend-item"><div class="legend-dot" style="background:#E84855"></div>Unassigned ZMWs</div>
  <div class="legend-item"><div class="legend-dot" style="background:#BBBBBB"></div>Dropout / Non-FLNC</div>
  <div class="legend-item"><div class="legend-dot" style="background:#457B9D"></div>FLNC polya reads</div>
  {''.join(f'<div class="legend-item"><div class="legend-dot" style="background:{c}"></div>{l}</div>' for l, c in lib_color_map.items())}
  <div style="width:100%;border-top:1px solid #e5e7eb;margin:0.25rem 0"></div>
  {''.join(f'<div class="legend-item"><div class="legend-dot" style="background:{c}"></div>{bc}</div>' for bc, c in BC_COLORS.items())}
</div>
<div id="sankey"></div>
<script>
  const fig = {fig_json};
  Plotly.newPlot('sankey', fig.data, fig.layout, {{
    responsive: true,
    displayModeBar: true,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'],
    toImageButtonOptions: {{format: 'png', filename: 'pipeline_read_flow', scale: 2}}
  }});
</script>
</body>
</html>"""

out_html.parent.mkdir(parents=True, exist_ok=True)
out_html.write_text(html)
print(f"Written to {out_html}")
