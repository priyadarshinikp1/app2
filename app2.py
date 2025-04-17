#python -m streamlit run app.py

"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from gseapy import enrichr

st.set_page_config(page_title="Multi-Omics App", layout="wide")
st.title("🧬 Multi-Omics Integration App")

# Upload files
genomics = st.file_uploader("Upload Genomics CSV", type="csv")
transcriptomics = st.file_uploader("Upload Transcriptomics CSV", type="csv")
proteomics = st.file_uploader("Upload Proteomics CSV", type="csv")

# Show previews
if genomics:
    gdf = pd.read_csv(genomics)
    st.subheader("Genomics Data")
    st.dataframe(gdf.head())

if transcriptomics:
    tdf = pd.read_csv(transcriptomics)
    st.subheader("Transcriptomics Data")
    st.dataframe(tdf.head())

if proteomics:
    pdf = pd.read_csv(proteomics)
    st.subheader("Proteomics Data")
    st.dataframe(pdf.head())

# ---- Filtering Section ----
st.header("🎛️ Filter & Integrate")

if genomics and transcriptomics and proteomics:
    st.subheader("🔬 Set Filtering Criteria")

    try:
        t_pval_thresh = float(st.text_input("Transcriptomics: Max p-value", value="0.05", key="t_pval"))
        p_intensity_thresh = float(st.text_input("Proteomics: Min Intensity", value="1000", key="p_intensity"))

        tdf_filtered = tdf[tdf['p_value'] <= t_pval_thresh]
        pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]
        gdf_filtered = gdf  # Genomics filter is paused

        st.write(f"✅ Transcriptomics filtered: {len(tdf_filtered)} genes")
        st.write(f"✅ Proteomics filtered: {len(pdf_filtered)} genes")
        st.write(f"✅ Genomics filtered (unchanged): {len(gdf_filtered)} genes")

        common_genes = set(gdf_filtered['Gene']) & set(tdf_filtered['Gene']) & set(pdf_filtered['Gene'])
        st.success(f"🎯 Genes common in all 3 omics: {len(common_genes)}")
        common_gene_list = list(common_genes)
        st.dataframe(pd.DataFrame({'Common Genes': common_gene_list}))

    except ValueError:
        st.error("❌ Please enter valid numeric values for filters.")

# ---- Enrichment Analysis Section ----
if 'common_gene_list' in locals() and len(common_gene_list) > 0:
    st.header("🧬 Enrichment Analysis with Enrichr")

    if st.button("Run Enrichment Analysis"):
        # PATHWAY ANALYSIS
        st.subheader("🧠 Enriched Pathways (Reactome) + Involved Genes")
        try:
            enr_path = enrichr(gene_list=common_gene_list,
                               gene_sets="Reactome_2016",
                               organism="Human",
                               outdir=None)

            df_pathways = enr_path.results
            df_pathways["-log10(pval)"] = -np.log10(df_pathways["P-value"])
            df_pathways = df_pathways.sort_values("P-value").reset_index(drop=True)

            # Table
            st.dataframe(df_pathways[['Term', 'P-value', 'Overlap', 'Genes']].rename(
                columns={"Term": "Pathway", "Genes": "Genes Involved"}))

            # Plot
            top_path = df_pathways.head(10)
            fig = px.bar(
                top_path,
                x="Term",
                y="-log10(pval)",
                title="Top 10 Enriched Pathways (Reactome)",
                labels={"Term": "Pathway", "-log10(pval)": "-log10(p)"},
            )
            fig.update_layout(xaxis_tickangle=-45)
            st.plotly_chart(fig, use_container_width=True)

        except Exception as e:
            st.error(f"Error in Reactome enrichment: {e}")

        # DISEASE ASSOCIATION
        st.subheader("🦠 Disease Associations (DisGeNET)")
        try:
            enr_disease = enrichr(gene_list=common_gene_list,
                                  gene_sets="DisGeNET",
                                  organism="Human",
                                  outdir=None)

            df_disease = enr_disease.results
            df_disease["-log10(pval)"] = -np.log10(df_disease["P-value"])
            df_disease = df_disease.sort_values("P-value").reset_index(drop=True)

            st.dataframe(df_disease[['Term', 'P-value', 'Genes']].rename(
                columns={"Term": "Disease", "Genes": "Associated Genes"}))

            top_disease = df_disease.head(10)
            fig = px.bar(
                top_disease,
                x="Term",
                y="-log10(pval)",
                title="Top 10 Associated Diseases",
                labels={"Term": "Disease", "-log10(pval)": "-log10(p)"},
            )
            fig.update_layout(xaxis_tickangle=-45)
            st.plotly_chart(fig, use_container_width=True)

        except Exception as e:
            st.error(f"Error in disease enrichment: {e}")

        # METABOLITE ASSOCIATION
        st.subheader("🧪 Metabolite Associations (HMDB Metabolites)")
        try:
            enr_metabo = enrichr(gene_list=common_gene_list,
                                 gene_sets="HMDB_Metabolites",
                                 organism="Human",
                                 outdir=None)

            df_metabo = enr_metabo.results
            df_metabo["-log10(pval)"] = -np.log10(df_metabo["P-value"])
            df_metabo = df_metabo.sort_values("P-value").reset_index(drop=True)

            st.dataframe(df_metabo[['Term', 'P-value', 'Genes']].rename(
                columns={"Term": "Metabolite", "Genes": "Associated Genes"}))

            top_metab = df_metabo.head(10)
            fig = px.bar(
                top_metab,
                x="Term",
                y="-log10(pval)",
                title="Top 10 Associated Metabolites",
                labels={"Term": "Metabolite", "-log10(pval)": "-log10(p)"},
            )
            fig.update_layout(xaxis_tickangle=-45)
            st.plotly_chart(fig, use_container_width=True)

        except Exception as e:
            st.error(f"Error in metabolite enrichment: {e}")
"""

"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from gseapy import enrichr

st.set_page_config(page_title="Multi-Omics App", layout="wide")
st.title("🧬 Multi-Omics Integration App")

# Upload files
genomics = st.file_uploader("Upload Genomics CSV", type="csv")
transcriptomics = st.file_uploader("Upload Transcriptomics CSV", type="csv")
proteomics = st.file_uploader("Upload Proteomics CSV", type="csv")

# Show previews
if genomics:
    gdf = pd.read_csv(genomics)
    st.subheader("Genomics Data")
    st.dataframe(gdf.head())

if transcriptomics:
    tdf = pd.read_csv(transcriptomics)
    st.subheader("Transcriptomics Data")
    st.dataframe(tdf.head())

if proteomics:
    pdf = pd.read_csv(proteomics)
    st.subheader("Proteomics Data")
    st.dataframe(pdf.head())

# ---- Filtering Section ----
st.header("🎛️ Filter & Integrate")

if genomics and transcriptomics and proteomics:
    st.subheader("🔬 Set Filtering Criteria")

    try:
        t_pval_thresh = float(st.text_input("Transcriptomics: Max p-value", value="0.05", key="t_pval"))
        p_intensity_thresh = float(st.text_input("Proteomics: Min Intensity", value="1000", key="p_intensity"))

        tdf_filtered = tdf[tdf['p_value'] <= t_pval_thresh]
        pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]
        gdf_filtered = gdf

        st.write(f"✅ Transcriptomics filtered: {len(tdf_filtered)} genes")
        st.write(f"✅ Proteomics filtered: {len(pdf_filtered)} genes")
        st.write(f"✅ Genomics filtered (unchanged): {len(gdf_filtered)} genes")

        common_genes = set(gdf_filtered['Gene']) & set(tdf_filtered['Gene']) & set(pdf_filtered['Gene'])
        st.success(f"🌯️️ Genes common in all 3 omics: {len(common_genes)}")
        common_genes_df = pd.DataFrame({'Common Genes': list(common_genes)})
        st.dataframe(common_genes_df)

    except ValueError:
        st.error("❌ Please enter valid numeric values for filters.")

    # Enrichment Analysis
    st.header("📈 Enrichment Analyses")
    if st.button("Run Enrichr Analyses"):
        libraries = {
            "Reactome Pathways": "Reactome_2016",
            "Disease Associations": "DisGeNET",
            "HMDB Metabolites": "HMDB_Metabolites"
        }

        results = {}
        edge_list = []
        assoc_data = []

        for name, lib in libraries.items():
            try:
                enr = enrichr(gene_list=list(common_genes), gene_sets=lib, outdir=None)
                if enr.results.empty:
                    st.warning(f"⚠️ No results from {name}")
                    continue
                df = enr.results.copy()
                df['-log10(pval)'] = -np.log10(df['P-value'])
                df = df.rename(columns={"Term": "Pathway", "Genes": "Genes_Involved"})
                results[name] = df

                st.subheader(f"📋 {name} - Top Results")
                st.dataframe(df[['Pathway', 'P-value', 'Adjusted P-value', 'Overlap', 'Genes_Involved']].head(10))

                fig = px.bar(
                    df.head(10),
                    x="Pathway",
                    y="-log10(pval)",
                    title=f"Top 10 {name}",
                    labels={"Pathway": "Term", "-log10(pval)": "-log10(p)"},
                )
                st.plotly_chart(fig)

            except Exception as e:
                st.error(f"Error in {name} enrichment: {e}")

        # 🔸 Network Visualization
        st.subheader("🔸 Omics-Disease-Metabolite Network")
        try:
            G = nx.Graph()
            color_map = {
                "Reactome Pathways": "lightblue",
                "Disease Associations": "lightcoral",
                "HMDB Metabolites": "lightgreen"
            }

            for name, df in results.items():
                color = color_map.get(name, "gray")
                for _, row in df.head(5).iterrows():
                    pathway = row['Pathway']
                    G.add_node(pathway, type=name, color=color)
                    for gene in row['Genes_Involved'].split(';'):
                        gene = gene.strip()
                        G.add_node(gene, type='Gene', color='lightgray')
                        G.add_edge(gene, pathway)
                        edge_list.append({'Gene': gene, 'Associated Term': pathway, 'Category': name})

                        # Prepare row for new association table
                        assoc_data.append({
                            'Gene': gene,
                            'Pathway': pathway if name == 'Reactome Pathways' else '',
                            'Metabolite': pathway if name == 'HMDB Metabolites' else '',
                            'Disease': pathway if name == 'Disease Associations' else ''
                        })

            pos = nx.spring_layout(G, seed=42)
            node_colors = [G.nodes[n]['color'] for n in G.nodes()]
            fig, ax = plt.subplots(figsize=(16, 12))
            nx.draw(
                G, pos,
                with_labels=True,
                node_color=node_colors,
                node_size=1200,
                font_size=10,
                font_weight='bold',
                edge_color='gray',
                ax=ax
            )

            # Legend
            legend_elements = [
                Patch(facecolor='lightblue', edgecolor='k', label='Reactome Pathways'),
                Patch(facecolor='lightcoral', edgecolor='k', label='Disease Associations'),
                Patch(facecolor='lightgreen', edgecolor='k', label='HMDB Metabolites'),
                Patch(facecolor='lightgray', edgecolor='k', label='Genes')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
            st.pyplot(fig)

            # Display association table
            edge_df = pd.DataFrame(edge_list)
            st.subheader("🔗 Gene-Term Associations Table")
            st.dataframe(edge_df)

            # Display organized association table
            assoc_df = pd.DataFrame(assoc_data)
            assoc_df = assoc_df.groupby(['Gene'], as_index=False).agg({
                'Pathway': lambda x: ', '.join(filter(None, x)),
                'Metabolite': lambda x: ', '.join(filter(None, x)),
                'Disease': lambda x: ', '.join(filter(None, x))
            })
            st.subheader("📄 Full Gene-Term Association Summary")
            st.dataframe(assoc_df)

        except Exception as e:
            st.error(f"Network drawing failed: {e}")
"""
"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from pyvis.network import Network
import networkx as nx
from gseapy import enrichr
import tempfile
import os
import requests

st.set_page_config(page_title="Multi-Omics App", layout="wide")
st.title("🧬 Multi-Omics Integration App")

# Upload files
genomics = st.file_uploader("Upload Genomics CSV", type="csv")
transcriptomics = st.file_uploader("Upload Transcriptomics CSV", type="csv")
proteomics = st.file_uploader("Upload Proteomics CSV", type="csv")

# Show previews
if genomics:
    gdf = pd.read_csv(genomics)
    st.subheader("Genomics Data")
    st.dataframe(gdf.head())

if transcriptomics:
    tdf = pd.read_csv(transcriptomics)
    st.subheader("Transcriptomics Data")
    st.dataframe(tdf.head())

if proteomics:
    pdf = pd.read_csv(proteomics)
    st.subheader("Proteomics Data")
    st.dataframe(pdf.head())

# ---- Sidebar Settings ----
st.sidebar.header("⚙️ Settings")

# Filter thresholds
st.sidebar.subheader("🔬 Filter Criteria")
cadd_thresh = float(st.sidebar.text_input("Min CADD Score (Genomics)", value="20"))
logfc_thresh = float(st.sidebar.text_input("Min |logFC| (Transcriptomics)", value="1"))
t_pval_thresh = float(st.sidebar.text_input("Max p-value (Transcriptomics)", value="0.05"))
p_intensity_thresh = float(st.sidebar.text_input("Min Intensity (Proteomics)", value="1000"))

# Display toggles
st.sidebar.subheader("🗽️ Display Controls")
run_enrichment = st.sidebar.checkbox("Run Enrichment Analyses", value=True)
show_network = st.sidebar.checkbox("Show Network Visualization", value=True)
show_association_table = st.sidebar.checkbox("Show Association Table", value=True)

# Pathway display control
num_pathways_to_show = int(st.sidebar.slider("Number of Pathways to Display in Network", min_value=1, max_value=50, value=10))

# ---- Filtering & Integration ----
st.header("🎛️ Filter & Integrate")

if genomics and transcriptomics and proteomics:
    try:
        if 'CADD' not in gdf.columns:
            st.error("❌ 'CADD' column not found in Genomics file.")
            st.stop()

        gdf_filtered = gdf[gdf['CADD'] >= cadd_thresh]
        tdf_filtered = tdf[(tdf['p_value'] <= t_pval_thresh) & (tdf['logFC'].abs() >= logfc_thresh)]
        pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]

        st.write(f"✅ Genomics filtered (CADD ≥ {cadd_thresh}): {len(gdf_filtered)} genes")
        st.write(f"✅ Transcriptomics filtered (|logFC| ≥ {logfc_thresh}, p ≤ {t_pval_thresh}): {len(tdf_filtered)} genes")
        st.write(f"✅ Proteomics filtered (Intensity ≥ {p_intensity_thresh}): {len(pdf_filtered)} proteins")

        union_genes = set(gdf_filtered['Gene']) | set(tdf_filtered['Gene'])
        proteins = set(pdf_filtered['Protein'])

        def map_uniprot_to_gene(uniprot_ids):
            mapping = {}
            ids = list(uniprot_ids)
            for i in range(0, len(ids), 100):
                chunk = ids[i:i+100]
                query = " OR ".join([f"accession:{id_}" for id_ in chunk])
                url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,gene_names&format=tsv"
                try:
                    r = requests.get(url)
                    if r.status_code == 200:
                        lines = r.text.strip().split('\n')[1:]
                        for line in lines:
                            acc, genes = line.split('\t')
                            mapping[acc] = genes.split()[0] if genes else acc
                except Exception as e:
                    print(f"Error with UniProt API chunk: {e}")
                    continue
            return mapping

        st.info("🔄 Mapping UniProt IDs to gene names via UniProt API...")
        uniprot_gene_map = map_uniprot_to_gene(proteins)
        pdf_filtered['GeneName'] = pdf_filtered['Protein'].map(uniprot_gene_map)
        protein_gene_map = dict(zip(pdf_filtered['Protein'], pdf_filtered['GeneName']))

        all_entities = union_genes | set(protein_gene_map.values())
        st.success(f"🔗 Genes/Proteins from any omics layer: {len(all_entities)}")
        st.dataframe(pd.DataFrame({'Genes/Proteins': list(all_entities)}))

        results = {}
        raw_assoc_data = []

        if run_enrichment:
            st.header("📊 Enrichment Analyses")
            libraries = {
                "Reactome Pathways": "Reactome_2016",
                "Disease Associations": "DisGeNET",
                "HMDB Metabolites": "HMDB_Metabolites"
            }

            for name, lib in libraries.items():
                try:
                    enr = enrichr(gene_list=list(union_genes), gene_sets=lib, outdir=None)
                    if enr.results.empty:
                        st.warning(f"⚠️ No results from {name}")
                        continue
                    df = enr.results.copy()
                    df['-log10(pval)'] = -np.log10(df['P-value'])
                    df = df.rename(columns={"Term": "Pathway", "Genes": "Genes_Involved"})
                    results[name] = df

                    st.subheader(f"📋 {name} - Top Results")
                    st.dataframe(df[['Pathway', 'P-value', 'Adjusted P-value', 'Overlap', 'Genes_Involved']].head(10))

                    fig = px.bar(
                        df.head(10),
                        x="Pathway",
                        y="-log10(pval)",
                        title=f"Top 10 {name}",
                        labels={"Pathway": "Term", "-log10(pval)": "-log10(p)"},
                    )
                    st.plotly_chart(fig)

                except Exception as e:
                    st.error(f"Error in {name} enrichment: {e}")

        if show_network and results:
            st.subheader("🧠 Interactive Omics Network")
            try:
                net = Network(height='800px', width='100%', notebook=False, directed=False)
                net.force_atlas_2based()
                net.show_buttons(filter_=['physics'])

                legend_items = {
                    "Gene": '#a9a9a9',
                    "Protein": '#ffd700',
                    "Pathway": '#87cefa',
                    "Metabolite": '#98fb98',
                    "Disease": '#f08080'
                }
                y_pos = 0
                for label, color in legend_items.items():
                    net.add_node(f"legend_{label}", label=label, shape='box', color=color, size=20, x=-1000, y=y_pos, physics=False, fixed=True)
                    y_pos -= 50

                color_map = {
                    "Reactome Pathways": "#87cefa",
                    "Disease Associations": "#f08080",
                    "HMDB Metabolites": "#98fb98"
                }

                for name, df in results.items():
                    color = color_map.get(name, "gray")
                    for _, row in df.head(num_pathways_to_show).iterrows():
                        term = row['Pathway']
                        net.add_node(term, label=term, color=color)
                        for gene in row['Genes_Involved'].split(';'):
                            gene = gene.strip()
                            if not gene:
                                continue
                            net.add_node(gene, label=gene, color='#a9a9a9')
                            net.add_edge(gene, term)

                            matched_proteins = [prot for prot, gname in protein_gene_map.items() if gname == gene]
                            for prot in matched_proteins:
                                net.add_node(prot, label=prot, color='#ffd700')
                                net.add_edge(gene, prot)

                            raw_assoc_data.append({
                                'Gene': gene,
                                'Protein': ';'.join(matched_proteins) if matched_proteins else '',
                                'Pathway': term if name == 'Reactome Pathways' else '',
                                'Metabolite': term if name == 'HMDB Metabolites' else '',
                                'Disease': term if name == 'Disease Associations' else ''
                            })

                with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
                    net.save_graph(tmp_file.name)
                    st.components.v1.html(open(tmp_file.name, 'r', encoding='utf-8').read(), height=800)

            except Exception as e:
                st.error(f"Network rendering failed: {e}")

        # ---- Association Table ----
        if show_association_table and raw_assoc_data:
            df = pd.DataFrame(raw_assoc_data)
            assoc_df = df.groupby('Gene').agg({
                'Protein': lambda x: ';'.join(set(filter(None, x))),
                'Pathway': lambda x: ';'.join(set(filter(None, x))),
                'Disease': lambda x: ';'.join(set(filter(None, x))),
                'Metabolite': lambda x: ';'.join(set(filter(None, x)))
            }).reset_index()

            assoc_df = assoc_df[['Disease', 'Pathway', 'Gene', 'Protein', 'Metabolite']]
            assoc_df['non_nulls'] = assoc_df.notnull().sum(axis=1)
            assoc_df = assoc_df.sort_values(by='non_nulls', ascending=False).drop(columns='non_nulls')

            st.subheader("📄 Gene-Protein-Term Association Summary")
            st.dataframe(assoc_df)

    except ValueError:
        st.error("❌ Please enter valid numeric values for filters.")
"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from pyvis.network import Network
import networkx as nx
from gseapy import enrichr
from sklearn.decomposition import PCA
import umap
import tempfile
import os
import requests

st.set_page_config(page_title="Multi-Omics App", layout="wide")
st.title("🧬 Multi-Omics Integration App")

# Upload files
genomics = st.file_uploader("Upload Genomics CSV", type="csv")
transcriptomics = st.file_uploader("Upload Transcriptomics CSV", type="csv")
proteomics = st.file_uploader("Upload Proteomics CSV", type="csv")

# Show previews
if genomics:
    gdf = pd.read_csv(genomics)
    st.subheader("Genomics Data")
    st.dataframe(gdf.head())

if transcriptomics:
    tdf = pd.read_csv(transcriptomics)
    st.subheader("Transcriptomics Data")
    st.dataframe(tdf.head())

if proteomics:
    pdf = pd.read_csv(proteomics)
    st.subheader("Proteomics Data")
    st.dataframe(pdf.head())

# Sidebar controls
st.sidebar.header("⚙️ Settings")
cadd_thresh = float(st.sidebar.text_input("Min CADD Score (Genomics)", value="20"))
logfc_thresh = float(st.sidebar.text_input("Min |logFC| (Transcriptomics)", value="1"))
t_pval_thresh = float(st.sidebar.text_input("Max p-value (Transcriptomics)", value="0.05"))
p_intensity_thresh = float(st.sidebar.text_input("Min Intensity (Proteomics)", value="1000"))

run_enrichment = st.sidebar.checkbox("Run Enrichment Analyses", value=True)
show_network = st.sidebar.checkbox("Show Network Visualization", value=True)
show_association_table = st.sidebar.checkbox("Show Association Table", value=True)

num_pathways_to_show = int(st.sidebar.slider("Number of Pathways to Display in Network", min_value=1, max_value=50, value=10))

# Filtering and Integration
st.header("🎛️ Filter & Integrate")

if genomics and transcriptomics and proteomics:
    try:
        gdf_filtered = gdf[gdf['CADD'] >= cadd_thresh]
        tdf_filtered = tdf[(tdf['p_value'] <= t_pval_thresh) & (tdf['logFC'].abs() >= logfc_thresh)]
        pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]

        st.write(f"✅ Genomics filtered: {len(gdf_filtered)}")
        st.write(f"✅ Transcriptomics filtered: {len(tdf_filtered)}")
        st.write(f"✅ Proteomics filtered: {len(pdf_filtered)}")

        union_genes = set(gdf_filtered['Gene']) | set(tdf_filtered['Gene'])
        proteins = set(pdf_filtered['Protein'])

        def map_uniprot_to_gene(uniprot_ids):
            mapping = {}
            ids = list(uniprot_ids)
            for i in range(0, len(ids), 100):
                chunk = ids[i:i+100]
                query = " OR ".join([f"accession:{id_}" for id_ in chunk])
                url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,gene_names&format=tsv"
                try:
                    r = requests.get(url)
                    if r.status_code == 200:
                        lines = r.text.strip().split('\n')[1:]
                        for line in lines:
                            acc, genes = line.split('\t')
                            mapping[acc] = genes.split()[0] if genes else acc
                except Exception as e:
                    print(f"Error with UniProt API chunk: {e}")
            return mapping

        st.info("🔄 Mapping UniProt IDs to gene names via UniProt API...")
        uniprot_gene_map = map_uniprot_to_gene(proteins)
        pdf_filtered['GeneName'] = pdf_filtered['Protein'].map(uniprot_gene_map)
        protein_gene_map = dict(zip(pdf_filtered['Protein'], pdf_filtered['GeneName']))

        all_entities = union_genes | set(protein_gene_map.values())
        st.success(f"🔗 Unique genes/proteins across layers: {len(all_entities)}")
        st.dataframe(pd.DataFrame({'Genes/Proteins': list(all_entities)}))

        results = {}
        raw_assoc_data = []

        if run_enrichment:
            st.header("📊 Enrichment Analyses")
            libraries = {
                "Reactome Pathways": "Reactome_2016",
                "Disease Associations": "DisGeNET",
                "HMDB Metabolites": "HMDB_Metabolites"
            }

            for name, lib in libraries.items():
                try:
                    gene_list_clean = [str(g).strip() for g in union_genes if pd.notna(g) and str(g).strip()]
                    enr = enrichr(gene_list=gene_list_clean, gene_sets=lib, outdir=None)
                    if enr.results.empty:
                        st.warning(f"⚠️ No results from {name}")
                        continue
                    df = enr.results.copy()
                    df['-log10(pval)'] = -np.log10(df['P-value'])
                    df = df.rename(columns={"Term": "Pathway", "Genes": "Genes_Involved"})
                    results[name] = df

                    st.subheader(f"📋 {name} - Top Results")
                    st.dataframe(df[['Pathway', 'P-value', 'Adjusted P-value', 'Overlap', 'Genes_Involved']].head(10))

                    fig = px.bar(
                        df.head(10),
                        x="Pathway",
                        y="-log10(pval)",
                        title=f"Top 10 {name}",
                        labels={"Pathway": "Term", "-log10(pval)": "-log10(p)"},
                    )
                    st.plotly_chart(fig)

                except Exception as e:
                    st.error(f"Error in {name} enrichment: {e}")

        if show_network and results:
            st.subheader("🧠 Interactive Omics Network")
            try:
                net = Network(height='800px', width='100%', notebook=False, directed=False)
                net.force_atlas_2based()
                net.show_buttons(filter_=['physics'])

                legend_items = {
                    "Gene": 'rgb(169,169,169)',
                    "Protein": 'rgb(255,215,0)',
                    "Pathway": 'rgb(135,206,250)',
                    "Metabolite": 'rgb(152,251,152)',
                    "Disease": 'rgb(240,128,128)'
                }

                y_pos = 0
                for label, color in legend_items.items():
                    net.add_node(f"legend_{label}", label=label, shape='box', color=color, size=20, x=-1000, y=y_pos, physics=False, fixed=True)
                    y_pos -= 50

                color_map = {
                    "Reactome Pathways": "rgb(135,206,250)",
                    "Disease Associations": "rgb(240,128,128)",
                    "HMDB Metabolites": "rgb(152,251,152)"
                }

                for name, df in results.items():
                    color = color_map.get(name, "gray")
                    for _, row in df.head(num_pathways_to_show).iterrows():
                        term = row['Pathway']
                        net.add_node(term, label=term, color=color)
                        for gene in row['Genes_Involved'].split(';'):
                            gene = gene.strip()
                            if not gene:
                                continue
                            net.add_node(gene, label=gene, color='rgb(169,169,169)')
                            net.add_edge(gene, term)

                            matched_proteins = [prot for prot, gname in protein_gene_map.items() if gname == gene]
                            for prot in matched_proteins:
                                net.add_node(prot, label=prot, color='rgb(255,215,0)')
                                net.add_edge(gene, prot)

                            raw_assoc_data.append({
                                'Gene': gene,
                                'Protein': ';'.join(matched_proteins) if matched_proteins else '',
                                'Pathway': term if name == 'Reactome Pathways' else '',
                                'Metabolite': term if name == 'HMDB Metabolites' else '',
                                'Disease': term if name == 'Disease Associations' else ''
                            })

                with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
                    net.save_graph(tmp_file.name)
                    st.components.v1.html(open(tmp_file.name, 'r', encoding='utf-8').read(), height=800)

            except Exception as e:
                st.error(f"Network rendering failed: {e}")

        if show_association_table and raw_assoc_data:
            df = pd.DataFrame(raw_assoc_data)
            assoc_df = df.groupby('Gene').agg({
                'Protein': lambda x: ';'.join(set(filter(None, x))),
                'Pathway': lambda x: ';'.join(set(filter(None, x))),
                'Disease': lambda x: ';'.join(set(filter(None, x))),
                'Metabolite': lambda x: ';'.join(set(filter(None, x)))
            }).reset_index()

            assoc_df['non_nulls'] = assoc_df.notnull().sum(axis=1)
            assoc_df = assoc_df.sort_values(by='non_nulls', ascending=False).drop(columns='non_nulls')

            st.subheader("📄 Gene-Protein-Term Association Summary")
            st.dataframe(assoc_df)

        # 🧪 Dimensionality Reduction: PCA + UMAP
        if transcriptomics:
            st.header("🔍 Transcriptomics: PCA and UMAP")
            try:
                expr_data = tdf.dropna()
                feature_cols = [col for col in expr_data.columns if expr_data[col].dtype != 'object']
                X = expr_data[feature_cols]

                pca = PCA(n_components=2)
                pca_result = pca.fit_transform(X)
                pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
                st.subheader("PCA Plot")
                st.plotly_chart(px.scatter(pca_df, x="PC1", y="PC2", title="PCA of Transcriptomics"))

                reducer = umap.UMAP(n_components=2)
                umap_result = reducer.fit_transform(X)
                umap_df = pd.DataFrame(umap_result, columns=["UMAP1", "UMAP2"])
                st.subheader("UMAP Plot")
                st.plotly_chart(px.scatter(umap_df, x="UMAP1", y="UMAP2", title="UMAP of Transcriptomics"))
            except Exception as e:
                st.warning(f"Error in dimensionality reduction: {e}")

    except ValueError:
        st.error("❌ Please enter valid numeric values for filters.")






















###########################################
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from pyvis.network import Network
import networkx as nx
from gseapy import enrichr
from sklearn.decomposition import PCA
import umap
import tempfile
import os
import requests

st.set_page_config(page_title="Multi-Omics App", layout="wide")
st.title("🧬 Multi-Omics Integration App")

# Upload files
genomics = st.file_uploader("Upload Genomics CSV", type="csv")
transcriptomics = st.file_uploader("Upload Transcriptomics CSV", type="csv")
proteomics = st.file_uploader("Upload Proteomics CSV", type="csv")

# Show previews
if genomics:
    gdf = pd.read_csv(genomics)
    st.subheader("Genomics Data")
    st.dataframe(gdf.head())

if transcriptomics:
    tdf = pd.read_csv(transcriptomics)
    st.subheader("Transcriptomics Data")
    st.dataframe(tdf.head())

if proteomics:
    pdf = pd.read_csv(proteomics)
    st.subheader("Proteomics Data")
    st.dataframe(pdf.head())

# Sidebar controls
st.sidebar.header("⚙️ Settings")
cadd_thresh = float(st.sidebar.text_input("Min CADD Score (Genomics)", value="20"))
logfc_thresh = float(st.sidebar.text_input("Min |logFC| (Transcriptomics)", value="1"))
t_pval_thresh = float(st.sidebar.text_input("Max p-value (Transcriptomics)", value="0.05"))
p_intensity_thresh = float(st.sidebar.text_input("Min Intensity (Proteomics)", value="1000"))

run_enrichment = st.sidebar.checkbox("Run Enrichment Analyses", value=True)
show_network = st.sidebar.checkbox("Show Network Visualization", value=True)
show_association_table = st.sidebar.checkbox("Show Association Table", value=True)

num_pathways_to_show = int(st.sidebar.slider("Number of Pathways to Display in Network", min_value=1, max_value=50, value=10))

# Filtering and Integration
st.header("🎛️ Filter & Integrate")

if genomics and transcriptomics and proteomics:
    try:
        gdf_filtered = gdf[gdf['CADD'] >= cadd_thresh]
        tdf_filtered = tdf[(tdf['p_value'] <= t_pval_thresh) & (tdf['logFC'].abs() >= logfc_thresh)]
        pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]

        st.write(f"✅ Genomics filtered: {len(gdf_filtered)}")
        st.write(f"✅ Transcriptomics filtered: {len(tdf_filtered)}")
        st.write(f"✅ Proteomics filtered: {len(pdf_filtered)}")

        union_genes = set(gdf_filtered['Gene']) | set(tdf_filtered['Gene'])

        def extract_uniprot_ids(protein_series):
            ids = set()
            for entry in protein_series.dropna():
                for pid in str(entry).split(";"):
                    pid = pid.strip()
                    if pid:
                        ids.add(pid)
            return ids

        def map_uniprot_to_gene(uniprot_ids):
            mapping = {}
            ids = list(uniprot_ids)
            for i in range(0, len(ids), 100):
                chunk = ids[i:i+100]
                query = " OR ".join([f"accession:{id_}" for id_ in chunk])
                url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,gene_names&format=tsv"
                try:
                    r = requests.get(url)
                    if r.status_code == 200:
                        lines = r.text.strip().split('\n')[1:]
                        for line in lines:
                            acc, genes = line.split('\t')
                            mapping[acc] = genes.split()[0] if genes else acc
                except Exception as e:
                    print(f"Error with UniProt API chunk: {e}")
            return mapping

        st.info("🔄 Mapping UniProt IDs to gene names via UniProt API...")
        unique_uniprot_ids = extract_uniprot_ids(pdf_filtered['Protein'])
        uniprot_gene_map = map_uniprot_to_gene(unique_uniprot_ids)

        expanded_rows = []
        for _, row in pdf_filtered.iterrows():
            protein_ids = str(row['Protein']).split(';')
            for pid in protein_ids:
                pid = pid.strip()
                if pid:
                    gene = uniprot_gene_map.get(pid, None)
                    if gene:
                        expanded_rows.append({'Protein': pid, 'GeneName': gene})

        expanded_protein_df = pd.DataFrame(expanded_rows)

        if not expanded_protein_df.empty and 'GeneName' in expanded_protein_df.columns:
            protein_gene_map = dict(zip(expanded_protein_df['Protein'], expanded_protein_df['GeneName']))
        else:
            st.warning("⚠️ No proteins could be mapped to gene names. Network may be incomplete.")
            protein_gene_map = {}

        st.write(f"🧪 Mapped {len(expanded_protein_df)} proteins to genes")
        st.dataframe(expanded_protein_df.head())

        all_entities = union_genes | set(protein_gene_map.values())
        st.success(f"🔗 Unique genes/proteins across layers: {len(all_entities)}")
        st.dataframe(pd.DataFrame({'Genes/Proteins': list(all_entities)}))

        results = {}
        raw_assoc_data = []

        if run_enrichment:
            st.header("📊 Enrichment Analyses")
            libraries = {
                "Reactome Pathways": "Reactome_2016",
                "Disease Associations": "DisGeNET",
                "HMDB Metabolites": "HMDB_Metabolites"
            }

            for name, lib in libraries.items():
                try:
                    gene_list_clean = [str(g).strip() for g in union_genes if pd.notna(g) and str(g).strip()]
                    enr = enrichr(gene_list=gene_list_clean, gene_sets=lib, outdir=None)
                    if enr.results.empty:
                        st.warning(f"⚠️ No results from {name}")
                        continue
                    df = enr.results.copy()
                    df['-log10(pval)'] = -np.log10(df['P-value'])
                    df = df.rename(columns={"Term": "Pathway", "Genes": "Genes_Involved"})
                    results[name] = df

                    st.subheader(f"📋 {name} - Top Results")
                    st.dataframe(df[['Pathway', 'P-value', 'Adjusted P-value', 'Overlap', 'Genes_Involved']].head(10))

                    fig = px.bar(
                        df.head(10),
                        x="Pathway",
                        y="-log10(pval)",
                        title=f"Top 10 {name}",
                        labels={"Pathway": "Term", "-log10(pval)": "-log10(p)"},
                    )
                    st.plotly_chart(fig)

                except Exception as e:
                    st.error(f"Error in {name} enrichment: {e}")

        if show_network and results:
            st.subheader("🧠 Interactive Omics Network")
            try:
                net = Network(height='800px', width='100%', notebook=False, directed=False)
                net.force_atlas_2based()
                net.show_buttons(filter_=['physics'])

                legend_items = {
                    "Gene": 'rgb(169,169,169)',
                    "Protein": 'rgb(255,215,0)',
                    "Pathway": 'rgb(135,206,250)',
                    "Metabolite": 'rgb(152,251,152)',
                    "Disease": 'rgb(240,128,128)'
                }

                y_pos = 0
                for label, color in legend_items.items():
                    net.add_node(f"legend_{label}", label=label, shape='box', color=color, size=20, x=-1000, y=y_pos, physics=False, fixed=True)
                    y_pos -= 50

                color_map = {
                    "Reactome Pathways": "rgb(135,206,250)",
                    "Disease Associations": "rgb(240,128,128)",
                    "HMDB Metabolites": "rgb(152,251,152)"
                }

                for name, df in results.items():
                    color = color_map.get(name, "gray")
                    for _, row in df.head(num_pathways_to_show).iterrows():
                        term = row['Pathway']
                        net.add_node(term, label=term, color=color)
                        for gene in row['Genes_Involved'].split(';'):
                            gene = gene.strip()
                            if not gene:
                                continue
                            net.add_node(gene, label=gene, color='rgb(169,169,169)')
                            net.add_edge(gene, term)

                            matched_proteins = [prot for prot, gname in protein_gene_map.items() if gname == gene]
                            for prot in matched_proteins:
                                net.add_node(prot, label=prot, color='rgb(255,215,0)')
                                net.add_edge(gene, prot)

                            raw_assoc_data.append({
                                'Gene': gene,
                                'Protein': ';'.join(matched_proteins) if matched_proteins else '',
                                'Pathway': term if name == 'Reactome Pathways' else '',
                                'Metabolite': term if name == 'HMDB Metabolites' else '',
                                'Disease': term if name == 'Disease Associations' else ''
                            })

                with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
                    net.save_graph(tmp_file.name)
                    st.components.v1.html(open(tmp_file.name, 'r', encoding='utf-8').read(), height=800)

            except Exception as e:
                st.error(f"Network rendering failed: {e}")

        if show_association_table and raw_assoc_data:
            df = pd.DataFrame(raw_assoc_data)
            assoc_df = df.groupby('Gene').agg({
                'Protein': lambda x: ';'.join(set(filter(None, x))),
                'Pathway': lambda x: ';'.join(set(filter(None, x))),
                'Disease': lambda x: ';'.join(set(filter(None, x))),
                'Metabolite': lambda x: ';'.join(set(filter(None, x)))
            }).reset_index()

            assoc_df['non_nulls'] = assoc_df.notnull().sum(axis=1)
            assoc_df = assoc_df.sort_values(by='non_nulls', ascending=False).drop(columns='non_nulls')

            st.subheader("📄 Gene-Protein-Term Association Summary")
            st.dataframe(assoc_df)

        # 🧪 Dimensionality Reduction: PCA + UMAP
        if transcriptomics:
            st.header("🔍 Transcriptomics: PCA and UMAP")
            try:
                expr_data = tdf.dropna()
                feature_cols = [col for col in expr_data.columns if expr_data[col].dtype != 'object']
                X = expr_data[feature_cols]

                if X.shape[0] > 2 and X.shape[1] > 2:
                    pca = PCA(n_components=2)
                    pca_result = pca.fit_transform(X)
                    pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
                    st.subheader("PCA Plot")
                    st.plotly_chart(px.scatter(pca_df, x="PC1", y="PC2", title="PCA of Transcriptomics"))

                    reducer = umap.UMAP(n_components=2)
                    umap_result = reducer.fit_transform(X)
                    umap_df = pd.DataFrame(umap_result, columns=["UMAP1", "UMAP2"])
                    st.subheader("UMAP Plot")
                    st.plotly_chart(px.scatter(umap_df, x="UMAP1", y="UMAP2", title="UMAP of Transcriptomics"))
                else:
                    st.warning("PCA/UMAP requires at least 3 samples and 3 features.")

            except Exception as e:
                st.warning(f"Error in dimensionality reduction: {e}")

    except ValueError:
        st.error("❌ Please enter valid numeric values for filters.")
