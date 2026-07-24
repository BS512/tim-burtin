import pandas as pd
import numpy as np
import altair as alt
import streamlit as st

# 1. Load Data
url = "https://raw.githubusercontent.com/BS512/tim-burtin/main/burtin.json"
df = pd.read_json(url)

# 2. Data Preprocessing
df_melted = df.melt(
    id_vars=["Bacteria", "Gram_Staining", "Genus"],
    value_vars=["Penicillin", "Streptomycin", "Neomycin"],
    var_name="Antibiotic",
    value_name="MIC"
)

df_melted["Gram_Staining"] = df_melted["Gram_Staining"].str.strip().str.lower()

# Logarithmic scaling & score calculation
df_melted["Log_MIC"] = np.log10(df_melted["MIC"])
min_log, max_log = df_melted["Log_MIC"].min(), df_melted["Log_MIC"].max()
df_melted["Effectiveness_Score"] = (max_log - df_melted["Log_MIC"]) / (max_log - min_log) * 100

# Sort species so Gram-negative and Gram-positive are grouped cleanly
df_melted = df_melted.sort_values(by=["Gram_Staining", "Bacteria"])
sorted_species = df_melted["Bacteria"].unique().tolist()

# 3. Pure Python Color Interpolation (Pink for Negative, Purple for Positive)
def hex_to_rgb(hex_str):
    hex_str = hex_str.lstrip('#')
    return tuple(int(hex_str[i:i+2], 16) for i in (0, 2, 4))

def rgb_to_hex(rgb):
    return '#{:02x}{:02x}{:02x}'.format(int(rgb[0]), int(rgb[1]), int(rgb[2]))

def get_cell_color(row):
    pct = max(0.0, min(1.0, row["Effectiveness_Score"] / 100.0))
    if row["Gram_Staining"] == "negative":
        # Pink Gradient: Light (#fce4ec) to Deep (#c2185b)
        c1, c2 = hex_to_rgb("#fce4ec"), hex_to_rgb("#c2185b")
    else:
        # Purple Gradient: Light (#f3e5f5) to Deep (#4a148c)
        c1, c2 = hex_to_rgb("#f3e5f5"), hex_to_rgb("#4a148c")
    
    interpolated = [c1[i] + (c2[i] - c1[i]) * pct for i in range(3)]
    return rgb_to_hex(interpolated)

df_melted["Color"] = df_melted.apply(get_cell_color, axis=1)

# 4. Build Altair Heatmap
heatmap = alt.Chart(df_melted).mark_rect(stroke='white', strokeWidth=1).encode(
    x=alt.X('Antibiotic:N', title=None, sort=['Penicillin', 'Streptomycin', 'Neomycin']),
    y=alt.Y('Bacteria:N', title=None, sort=sorted_species),
    color=alt.Color('Color:N', scale=None), # Uses computed hex values directly
    tooltip=[
        alt.Tooltip('Bacteria:N', title='Bacteria'),
        alt.Tooltip('Gram_Staining:N', title='Gram Staining'),
        alt.Tooltip('Antibiotic:N', title='Antibiotic'),
        alt.Tooltip('MIC:Q', title='MIC (µg/ml)'),
        alt.Tooltip('Effectiveness_Score:Q', title='Effectiveness (%)', format='.1f')
    ]
).properties(
    width=380,
    height=520,
    title=alt.TitleParams(
        text="Antibiotic Effectiveness Matrix by Gram Staining",
        subtitle="Gram-Negative (Pink Gradient) vs Gram-Positive (Purple Gradient)",
        anchor="start",
        fontSize=16,
        subtitleFontSize=12
    )
)

# Rule separator dividing Gram groups
neg_count = df[df["Gram_Staining"].str.strip().str.lower() == "negative"]["Bacteria"].nunique()

divider_df = pd.DataFrame({'y': [sorted_species[neg_count - 1]]})
divider = alt.Chart(divider_df).mark_rule(
    color='gray', strokeWidth=2, strokeDash=[4, 4]
).encode(y=alt.Y('y:N', sort=sorted_species))

# Final Layout Assembly
final_chart = alt.layer(heatmap, divider).configure_view(strokeWidth=0)

# 5. Streamlit App Display & Custom Legend
st.set_page_config(page_title="Burtin Antibiotic Analysis", layout="centered")
st.title("Antibiotic Resistance & Gram Staining Correlation")

st.altair_chart(final_chart, use_container_width=True)

# Custom Streamlit HTML Legend (Clean & Minimalist)
st.markdown("""
<div style="display: flex; gap: 30px; font-size: 14px; margin-top: -10px;">
    <div>
        <span style="font-weight: bold; color: #c2185b;">■ Gram-Negative (Pink):</span> 
        Light Pink (Low Effectiveness) → Dark Pink (High Effectiveness)
    </div>
    <div>
        <span style="font-weight: bold; color: #4a148c;">■ Gram-Positive (Purple):</span> 
        Light Purple (Low Effectiveness) → Dark Purple (High Effectiveness)
    </div>
</div>
""", unsafe_allow_html=True)
