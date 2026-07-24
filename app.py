import pandas as pd
import numpy as np
import altair as alt
import streamlit as st

# 1. Load Data directly from the JSON source
url = "https://raw.githubusercontent.com/BS512/tim-burtin/main/burtin.json"
df = pd.read_json(url)

# 2. Data Preprocessing & Normalization
# Melt antibiotics into long format
df_melted = df.melt(
    id_vars=["Bacteria", "Gram_Staining", "Genus"],
    value_vars=["Penicillin", "Streptomycin", "Neomycin"],
    var_name="Antibiotic",
    value_name="MIC"
)

# Clean up Gram_Staining column strings
df_melted["Gram_Staining"] = df_melted["Gram_Staining"].str.strip().str.lower()

# Invert logarithmic MIC values so higher values represent higher effectiveness
df_melted["Log_MIC"] = np.log10(df_melted["MIC"])
min_log, max_log = df_melted["Log_MIC"].min(), df_melted["Log_MIC"].max()

# Percentage score (100% = most effective / lowest MIC; 0% = least effective / highest MIC)
df_melted["Effectiveness_Score"] = (max_log - df_melted["Log_MIC"]) / (max_log - min_log) * 100

# Group species by Gram Staining and sort alphabetically
df_melted = df_melted.sort_values(by=["Gram_Staining", "Bacteria"])
sorted_species = df_melted["Bacteria"].unique().tolist()

# Gram-negative species count to calculate separation line position
neg_species_count = df[df["Gram_Staining"].str.strip().str.lower() == "negative"]["Bacteria"].nunique()

# 3. Build Heatmap Layers with Independent Color Scales

# Base encoding shared across both layers
base = alt.Chart(df_melted).encode(
    x=alt.X('Antibiotic:N', title=None, sort=['Penicillin', 'Streptomycin', 'Neomycin']),
    y=alt.Y('Bacteria:N', title=None, sort=sorted_species),
    tooltip=[
        alt.Tooltip('Bacteria:N', title='Bacteria Species'),
        alt.Tooltip('Gram_Staining:N', title='Gram Staining'),
        alt.Tooltip('Antibiotic:N', title='Antibiotic'),
        alt.Tooltip('MIC:Q', title='MIC (µg/ml)'),
        alt.Tooltip('Effectiveness_Score:Q', title='Effectiveness Score (%)', format='.1f')
    ]
)

# Layer 1: Gram-Negative -> PINK palette
heatmap_neg = base.transform_filter(
    alt.datum.Gram_Staining == 'negative'
).mark_rect().encode(
    color=alt.Color(
        'Effectiveness_Score:Q',
        title="Gram-Neg Effectiveness (%)",
        scale=alt.Scale(
            domain=[0, 100],
            range=['#fce4ec', '#c2185b'] # Light pink to Deep Pink
        )
    )
)

# Layer 2: Gram-Positive -> PURPLE palette
heatmap_pos = base.transform_filter(
    alt.datum.Gram_Staining == 'positive'
).mark_rect().encode(
    color=alt.Color(
        'Effectiveness_Score:Q',
        title="Gram-Pos Effectiveness (%)",
        scale=alt.Scale(
            domain=[0, 100],
            range=['#f3e5f5', '#4a148c'] # Light purple to Deep Purple
        )
    )
)

# Rule separator line dividing Gram-negative and Gram-positive species
divider_line = alt.Chart(pd.DataFrame({'y': [neg_species_count - 0.5]})).mark_rule(
    color='black',
    strokeWidth=2,
    strokeDash=[4, 4]
).encode(
    y='y:Q'
)

# Text annotations
annotation_neg = alt.Chart(pd.DataFrame({'text': ['Gram-Negative']})).mark_text(
    align='left', baseline='top', dx=130, dy=-240, fontWeight='bold', color='#c2185b' # Pink
).encode(text='text:N')

annotation_pos = alt.Chart(pd.DataFrame({'text': ['Gram-Positive']})).mark_text(
    align='left', baseline='bottom', dx=130, dy=230, fontWeight='bold', color='#4a148c' # Purple
).encode(text='text:N')

# Layer all components and resolve color scales independently
final_chart = alt.layer(
    heatmap_neg, heatmap_pos, divider_line, annotation_neg, annotation_pos
).resolve_scale(
    color='independent'
).properties(
    width=350,
    height=550,
    title=alt.TitleParams(
        text="Antibiotic Effectiveness Matrix by Gram Staining",
        subtitle="Gram-Negative (Pink) vs Gram-Positive (Purple)",
        anchor="start",
        fontSize=18,
        subtitleFontSize=13
    )
).configure_view(
    strokeWidth=0
).configure_axis(
    labelFontSize=11,
    titleFontSize=12
)

# 4. Streamlit App Layout
st.set_page_config(page_title="Burtin Antibiotic Analysis", layout="centered")
st.title("Antibiotic Resistance & Gram Staining Correlation")
st.altair_chart(final_chart, use_container_width=True)
