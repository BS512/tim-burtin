import pandas as pd
import numpy as np
import altair as alt
import streamlit as st

# 1. Load Data directly from the JSON source or local dictionary
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

# Invert logarithmic MIC values so higher values represent higher effectiveness (suppression power)
# Using Log10 transform since MIC values range across several orders of magnitude (0.001 to 870)
df_melted["Log_MIC"] = np.log10(df_melted["MIC"])
min_log, max_log = df_melted["Log_MIC"].min(), df_melted["Log_MIC"].max()

# Percentage score relative to worst/best effectiveness across the dataset
# (100% = most effective / lowest MIC; 0% = least effective / highest MIC)
df_melted["Effectiveness_Score"] = (max_log - df_melted["Log_MIC"]) / (max_log - min_log) * 100

# Group species by Gram Staining to ensure clear separation
# Order bacteria alphabetically within each Gram status group
df_melted["Gram_Staining"] = df_melted["Gram_Staining"].str.strip().str.lower()
df_melted = df_melted.sort_values(by=["Gram_Staining", "Bacteria"])

# Gram-negative species count to calculate separation line position
neg_species_count = df[df["Gram_Staining"].str.strip().str.lower() == "negative"]["Bacteria"].nunique()

# 3. Build Altair Visualizations
# Heatmap base tile chart
heatmap = alt.Chart(df_melted).mark_rect().encode(
    x=alt.X('Antibiotic:N', title=None, sort=['Penicillin', 'Streptomycin', 'Neomycin']),
    y=alt.Y('Bacteria:N', title=None, sort=df_melted['Bacteria'].unique().tolist()),
    color=alt.Color(
        'Effectiveness_Score:Q',
        title="Effectiveness (%)",
        scale=alt.Scale(
            domain=[0, 100],
            # Mapping distinct color ramps for storytelling: Purple hue for Gram-Negative, Pink hue for Gram-Positive
            range=['#f3e5f5', '#4a148c'] # Blended purple palette
        )
    ),
    tooltip=[
        alt.Tooltip('Bacteria:N', title='Bacteria Species'),
        alt.Tooltip('Gram_Staining:N', title='Gram Staining'),
        alt.Tooltip('Antibiotic:N', title='Antibiotic'),
        alt.Tooltip('MIC:Q', title='MIC (µg/ml)'),
        alt.Tooltip('Effectiveness_Score:Q', title='Effectiveness Score (%)', format='.1f')
    ]
).properties(
    width=350,
    height=550,
    title=alt.TitleParams(
        text="Antibiotic Effectiveness Matrix by Gram Staining",
        subtitle="Comparing Minimum Inhibitory Concentration (MIC) Effectiveness across Bacteria",
        anchor="start",
        fontSize=18,
        subtitleFontSize=13
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

# Text annotation overlaying the chart
annotation_neg = alt.Chart(pd.DataFrame({'text': ['Gram-Negative Group']})).mark_text(
    align='left', baseline='top', dx=130, dy=-240, fontWeight='bold', color='#4a148c'
).encode(text='text:N')

annotation_pos = alt.Chart(pd.DataFrame({'text': ['Gram-Positive Group']})).mark_text(
    align='left', baseline='bottom', dx=130, dy=230, fontWeight='bold', color='#c2185b'
).encode(text='text:N')

# Layer components together
final_chart = alt.layer(heatmap, divider_line, annotation_neg, annotation_pos).configure_view(
    strokeWidth=0
).configure_axis(
    labelFontSize=11,
    titleFontSize=12
)

# 4. Streamlit App Layout
st.set_page_config(page_title="Burtin Antibiotic Analysis", layout="centered")
st.title("Antibiotic Resistance & Gram Staining Correlation")
st.altair_chart(final_chart, use_container_width=True)
