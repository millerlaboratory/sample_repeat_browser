from pathlib import Path
import pandas as pd
import seaborn as sns
from shiny import reactive, render
from shiny.express import input, ui
import matplotlib.pyplot as plt
import plotly.express as px
from shinywidgets import render_widget

sns.set_theme(style="white")
df = pd.read_csv(Path(__file__).parent / "data" /"first_100_str_data_count.tsv", sep='\t')
df2= pd.read_csv(Path(__file__).parent / "data" /"first_100_str_data_motif.tsv", sep='\t')
str = pd.DataFrame(df)
motif = pd.DataFrame(df2)

# Extract unique sample options from the 'Sample' column
id_options = df['disease'].unique()

# Create a dictionary for the select input
disease_dict = {disease: disease for disease in id_options}

ui.page_opts(title="1KGP ONT Tandem Repeat Browser", fillable=True)

with ui.sidebar():
    ui.input_select(
        "select",  
        "Select a disease:",
        disease_dict
    )
    ui.input_slider("bin_number", "Binwidth", min=1, max=50, value=1, step=1)

    "This app was built from data made avaible by the 1000 Genomes ONT Sequencing Consortium. Genontype data was generated using Vamos (Ren et al., 2023) and repeat locus metadata was sourced from the STRChive (https://strchive.org/index)"


with ui.layout_columns(fill=False):

    with ui.value_box():
        "Gene"
        @render.express
        def select_gene():
            gene_data = filtered_data()
            first_gene = gene_data['gene'].iloc[0]
            first_gene
    
    with ui.value_box():
        "Type"
        @render.express
        def select_type():
            type_data = filtered_data()
            first_type = type_data['type'].iloc[0]
            first_type

    with ui.value_box():
        "Locus Structure"
        @render.express
        def select_locus():
            locus_data = filtered_data()
            first_locus = locus_data['locus_structure'].iloc[0]
            first_locus
    
    with ui.value_box():
        "Inheritance"
        @render.express
        def select_inheritance():
            in_data = filtered_data()
            first_value = in_data['Inheritance'].iloc[0]
            first_value
    
    with ui.value_box():
        "Pathogenic Min"
        @render.express
        def select_min():
            min_data = filtered_data()
            mean_value = min_data['pathogenic_min'].mean()
            mean_value

    with ui.value_box():
        "Pathogenic Max"
        @render.express
        def select_max():
            max_data = filtered_data()
            mean_value_max = max_data['pathogenic_max'].mean()
            mean_value_max

with ui.layout_columns():
    with ui.card(full_screen=True):
        with ui.card_header(class_="d-flex justify-content-between align-items-center"):
            "Repeat Size Distribution"

        @render_widget
        def histogram():
            d=filtered_data()
            # Define the bin width
            bin_width = input.bin_number()
            fighist = px.histogram(d, x="count", color_discrete_sequence=['green'], nbins=int((d['count'].max() - d['count'].min()) / bin_width))
            fighist.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            showlegend=False
            )

            # Set x and y axis titles
            fighist.update_xaxes(title_text='Motif count')
            fighist.update_yaxes(title_text='Number of alleles')

            return fighist
    

    with ui.card(full_screen=True):
        with ui.card_header(class_="d-flex justify-content-between align-items-center"):
            "Allele Sequences"
    
        @render.plot
        def waterfall():
            wf=filtered_motif()
            wf['motif'] = wf['motif'].astype('category')
            wf['motif_codes'] = wf['motif'].cat.codes
            lengths = wf.groupby('sample_allele')['anno'].count()
            sorted_sample_alleles = lengths.sort_values(ascending=False).index
            wf['sample_allele'] = pd.Categorical(wf['sample_allele'], categories=sorted_sample_alleles, ordered=True)
            pivot_df = wf.pivot(index='sample_allele', columns='pos', values='motif_codes')
            fig2, ax = plt.subplots(figsize=(10, 10))
            sns.heatmap(pivot_df, cmap='viridis', ax=ax, cbar=False)

            # Remove x-axis labels and change y-axis label
            ax.set(ylabel='Alleles', xlabel='Position')
            ax.set_yticklabels([])

            # Create a color map normalized to the range of motif_codes
            cmap = plt.cm.viridis
            norm = plt.Normalize(wf['motif_codes'].min(), wf['motif_codes'].max())

            # Create a legend with custom patches
            import matplotlib.patches as mpatches
            motifs = wf['motif'].unique()
            motif_codes = wf['motif_codes'].unique()
            patches = [mpatches.Patch(color=cmap(norm(code)), label=motif) for code, motif in zip(motif_codes, motifs)]
            plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

            return fig2
    
        
    with ui.card(full_screen=True):
        ui.card_header("Data")
        @render.data_frame  
        def table():
            data = filtered_data()
            # Select the columns you want to display
            subset_data = data[['sample', 'count', 'length', 'allele', 'sex']]
            return render.DataGrid(subset_data, filters=True)

@reactive.calc
def filtered_data():
    selected_id = input.select()
    return str[str['disease'] == selected_id]

@reactive.calc
def filtered_motif():
    selected_id = input.select()
    return motif[motif['disease'] == selected_id]