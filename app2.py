import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import io
import requests

# Helper to load CSV with encoding handling and clean column names
def load_csv(file):
    try:
        df = pd.read_csv(io.StringIO(file.getvalue().decode('utf-8')))
    except UnicodeDecodeError:
        df = pd.read_csv(io.StringIO(file.getvalue().decode('latin1')))
    df.columns = df.columns.str.strip()  # Clean column names (remove extra spaces)
    return df

# Function to fetch metabolites from HMDB based on gene
def fetch_metabolites_from_hmdb(gene):
    url = f"https://hmdb.ca/api/1/metabolites/gene/{gene}"  # This is a placeholder URL; update with correct endpoint if available
    response = requests.get(url)
    if response.status_code == 200:
        metabolites = response.json()  # Assuming the response is in JSON format
        return metabolites
    else:
        return []

# Function to fetch pathways and diseases from SMPDB based on metabolite
def fetch_pathways_diseases_from_smpdb(metabolite):
    url = f"https://smpdb.ca/api/pathways/metabolite/{metabolite}"  # This is a placeholder URL; update with correct endpoint if available
    response = requests.get(url)
    if response.status_code == 200:
        pathways_diseases = response.json()  # Assuming the response is in JSON format
        return pathways_diseases
    else:
        return []

# Sidebar for file uploads
st.sidebar.header("Upload your CSV files")
genomics_file = st.sidebar.file_uploader("Upload genomics CSV", type=["csv"])
transcriptomics_file = st.sidebar.file_uploader("Upload transcriptomics CSV", type=["csv"])
proteomics_file = st.sidebar.file_uploader("Upload proteomics CSV", type=["csv"])
gene_association_file = st.sidebar.file_uploader("Upload gene associations CSV", type=["csv"])
abberant_enzyme_file = st.sidebar.file_uploader("Upload abberant enzyme associations CSV", type=["csv"])

if genomics_file and transcriptomics_file and proteomics_file:
    # Load files
    genomics_df = load_csv(genomics_file)
    transcriptomics_df = load_csv(transcriptomics_file)
    proteomics_df = load_csv(proteomics_file)
    gene_assoc_df = load_csv(gene_association_file)
    abberant_enzyme_df = load_csv(abberant_enzyme_file)
    
    # Extract common genes
    common_genes = set(genomics_df['Gene']).intersection(set(transcriptomics_df['Gene'])).intersection(set(proteomics_df['Gene']))
    
    # Gene-enzyme associations from gene_association_file and abberant enzyme data
    gene_links = gene_assoc_df[gene_assoc_df['Gene'].isin(common_genes)]
    abberant_enzymes = abberant_enzyme_df['Enzyme'].dropna().unique()

    # Create graph
    G = nx.Graph()

    # Add nodes for genes and enzymes
    for gene in common_genes:
        G.add_node(gene, type='gene')  # Adding 'type' as 'gene' for each common gene
    for enzyme in abberant_enzymes:
        G.add_node(enzyme, type='enzyme')  # Adding 'type' as 'enzyme' for each abberant enzyme

    # Add edges (gene-to-enzyme association)
    for _, row in gene_links.iterrows():
        G.add_edge(row['Gene'], row['Enzyme'], type='gene-enzyme')
    
    # Fetch metabolites for each common gene from HMDB
    for gene in common_genes:
        metabolites = fetch_metabolites_from_hmdb(gene)
        for metabolite in metabolites:
            G.add_node(metabolite, type='metabolite')
            G.add_edge(gene, metabolite, type='gene-metabolite')

            # Fetch associated diseases and pathways from SMPDB
            pathways_diseases = fetch_pathways_diseases_from_smpdb(metabolite)
            for pathway_disease in pathways_diseases:
                pathway = pathway_disease.get('Pathway')
                disease = pathway_disease.get('Disease')
                if pathway:
                    if pathway not in G:
                        G.add_node(pathway, type='pathway')
                    G.add_edge(metabolite, pathway, type='metabolite-pathway')
                if disease:
                    if disease not in G:
                        G.add_node(disease, type='disease')
                    G.add_edge(metabolite, disease, type='metabolite-disease')

    # Add enzyme-disease associations (from gene_association_file or abberant_enzyme_file)
    for _, row in gene_links.iterrows():
        G.add_edge(row['Enzyme'], row['Disease'], type='enzyme-disease')

    # Define colors for node types
    type_colors = {
        'gene': 'blue',
        'enzyme': 'green',
        'metabolite': 'orange',
        'pathway': 'purple',
        'disease': 'red'
    }

    # Create a plotly network plot
    node_x = []
    node_y = []
    node_color = []
    node_size = []

    # Position nodes using networkx spring layout
    pos = nx.spring_layout(G, seed=42)

    # Collect node details for visualization
    for node, (x, y) in pos.items():
        node_x.append(x)
        node_y.append(y)
        node_type = G.nodes[node].get('type', 'unknown')
        node_color.append(type_colors.get(node_type, 'gray'))
        node_size.append(10)

    # Define the edge list for Plotly visualization
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_y.append(y0)
        edge_y.append(y1)

    # Create scatter plot for nodes
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        marker=dict(
            showscale=True,
            colorscale='Viridis',
            size=node_size,
            color=node_color,
            line_width=2
        ),
        text=list(G.nodes()),  # Tooltip text for nodes
        hoverinfo='text'
    )

    # Create line plot for edges
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        mode='lines',
        line=dict(width=0.5, color='#888'),
        hoverinfo='none'
    )

    # Create figure
    fig = go.Figure(data=[edge_trace, node_trace])

    # Network Visualization Layout
    fig.update_layout(
        title=dict(text="Gene-Enzyme-Metabolite-Pathway-Disease Network", font=dict(size=18)),
        showlegend=True,
        hovermode='closest',
        margin=dict(b=20, l=5, r=5, t=40),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
    )

    # Display the network graph
    st.plotly_chart(fig)

    # Document interactions in a table
    interactions = []
    for edge in G.edges(data=True):
        interactions.append({
            'Node 1': edge[0],
            'Node 2': edge[1],
            'Interaction Type': edge[2].get('type', 'unknown')
        })
    interactions_df = pd.DataFrame(interactions)
    st.subheader("Interactions Table")
    st.dataframe(interactions_df)

else:
    st.sidebar.warning("Please upload all necessary files to generate the network.")
