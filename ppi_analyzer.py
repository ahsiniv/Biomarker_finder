"""
Module 4: Protein-Protein Interaction Network Analysis
Uses STRING database API to build PPI networks and find hub genes.
"""

import requests
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import json
import time
import warnings
warnings.filterwarnings('ignore')


class PPINetworkAnalyzer:
    def __init__(self, output_dir: str = "./outputs"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.string_api = "https://string-db.org/api"
        self.species = 9606  # Human
        
    def get_string_interactions(self, gene_list: list, 
                                 min_score: int = 400) -> pd.DataFrame:
        """
        Fetch PPI interactions from STRING database.
        
        Args:
            gene_list: List of gene symbols
            min_score: Minimum interaction score (0-1000)
            
        Returns:
            DataFrame of interactions
        """
        print(f"\n🕸️  Building PPI network with STRING database...")
        print(f"  Querying {len(gene_list)} genes (min score: {min_score})...")
        
        if not gene_list:
            return pd.DataFrame()
        
        # Limit to first 200 genes to avoid API limits
        query_genes = gene_list[:200]
        genes_str = "%0d".join(query_genes)
        
        try:
            url = f"{self.string_api}/json/network"
            params = {
                "identifiers": genes_str,
                "species": self.species,
                "required_score": min_score,
                "caller_identity": "biomarker_finder"
            }
            
            response = requests.post(url, data=params, timeout=60)
            
            if response.status_code != 200:
                print(f"  ✗ STRING API error: {response.status_code}")
                return pd.DataFrame()
            
            interactions = response.json()
            
            if not interactions:
                print("  ⚠️  No interactions found, trying lower threshold...")
                params["required_score"] = 200
                response = requests.post(url, data=params, timeout=60)
                interactions = response.json()
            
            if not interactions:
                return pd.DataFrame()
            
            # Parse interactions
            rows = []
            for item in interactions:
                rows.append({
                    "gene_a": item.get("preferredName_A", ""),
                    "gene_b": item.get("preferredName_B", ""),
                    "score": item.get("score", 0),
                    "nscore": item.get("nscore", 0),
                    "fscore": item.get("fscore", 0),
                    "escore": item.get("escore", 0),
                    "ascore": item.get("ascore", 0)
                })
            
            df = pd.DataFrame(rows)
            df = df[df["gene_a"] != df["gene_b"]]  # Remove self-loops
            
            print(f"  ✓ Found {len(df)} interactions between {len(gene_list)} genes")
            return df
            
        except Exception as e:
            print(f"  ✗ STRING query error: {e}")
            return pd.DataFrame()

    def build_network(self, interactions_df: pd.DataFrame) -> nx.Graph:
        """Build NetworkX graph from interactions."""
        G = nx.Graph()
        
        if interactions_df.empty:
            return G
        
        for _, row in interactions_df.iterrows():
            G.add_edge(
                row["gene_a"], 
                row["gene_b"],
                weight=row["score"]
            )
        
        return G

    def calculate_hub_genes(self, G: nx.Graph, top_n: int = 20) -> pd.DataFrame:
        """
        Calculate centrality metrics and identify hub genes.
        
        Args:
            G: NetworkX graph
            top_n: Number of top hub genes to return
            
        Returns:
            DataFrame with hub gene scores
        """
        if G.number_of_nodes() == 0:
            return pd.DataFrame()
        
        print(f"  📐 Calculating centrality for {G.number_of_nodes()} nodes...")
        
        metrics = {}
        
        # Degree centrality
        degree_centrality = nx.degree_centrality(G)
        degree_count = dict(G.degree())
        
        # Betweenness centrality (skip for large graphs)
        if G.number_of_nodes() <= 500:
            betweenness = nx.betweenness_centrality(G, normalized=True)
        else:
            betweenness = {n: 0 for n in G.nodes()}
        
        # Closeness centrality
        try:
            closeness = nx.closeness_centrality(G)
        except Exception:
            closeness = {n: 0 for n in G.nodes()}
        
        # PageRank (good for finding important nodes)
        try:
            pagerank = nx.pagerank(G, alpha=0.85)
        except Exception:
            pagerank = {n: 1/G.number_of_nodes() for n in G.nodes()}
        
        # Combined hub score (MCC-like)
        rows = []
        for gene in G.nodes():
            degree = degree_count.get(gene, 0)
            
            # Maximal Clique Centrality approximation
            try:
                ego_graph = nx.ego_graph(G, gene, radius=1)
                subgraph_edges = ego_graph.number_of_edges()
                mcc_score = subgraph_edges
            except Exception:
                mcc_score = degree
            
            rows.append({
                "gene": gene,
                "degree": degree,
                "betweenness": round(betweenness.get(gene, 0), 4),
                "closeness": round(closeness.get(gene, 0), 4),
                "pagerank": round(pagerank.get(gene, 0), 6),
                "mcc": mcc_score
            })
        
        hub_df = pd.DataFrame(rows)
        
        # Normalize and compute composite score
        for col in ["degree", "betweenness", "closeness", "pagerank", "mcc"]:
            max_val = hub_df[col].max()
            if max_val > 0:
                hub_df[f"{col}_norm"] = hub_df[col] / max_val
            else:
                hub_df[f"{col}_norm"] = 0
        
        hub_df["hub_score"] = (
            hub_df["degree_norm"] * 0.3 +
            hub_df["betweenness_norm"] * 0.2 +
            hub_df["closeness_norm"] * 0.15 +
            hub_df["pagerank_norm"] * 0.2 +
            hub_df["mcc_norm"] * 0.15
        )
        
        hub_df = hub_df.sort_values("hub_score", ascending=False)
        
        top_hubs = hub_df.head(top_n)
        print(f"  🏆 Top 5 hub genes: {', '.join(top_hubs.head(5)['gene'].tolist())}")
        
        return top_hubs

    def visualize_network(self, G: nx.Graph, hub_genes: pd.DataFrame, 
                          disease_name: str, top_n: int = 50):
        """Generate PPI network visualization."""
        if G.number_of_nodes() == 0:
            return
        
        try:
            # Use subgraph of top hub genes for visualization
            top_genes = hub_genes.head(top_n)["gene"].tolist() if not hub_genes.empty else list(G.nodes())[:top_n]
            subgraph = G.subgraph([n for n in top_genes if n in G.nodes()])
            
            if subgraph.number_of_nodes() == 0:
                subgraph = G.subgraph(list(G.nodes())[:top_n])
            
            fig, ax = plt.subplots(figsize=(16, 14))
            
            # Layout
            pos = nx.spring_layout(subgraph, k=2.5, seed=42, iterations=100)
            
            # Node sizes based on degree
            degrees = dict(subgraph.degree())
            max_deg = max(degrees.values()) if degrees else 1
            node_sizes = [300 + (degrees.get(n, 0) / max_deg) * 1500 for n in subgraph.nodes()]
            
            # Node colors based on hub score
            if not hub_genes.empty:
                hub_dict = dict(zip(hub_genes["gene"], hub_genes["hub_score"]))
                node_colors = [hub_dict.get(n, 0.1) for n in subgraph.nodes()]
            else:
                node_colors = [degrees.get(n, 0) for n in subgraph.nodes()]
            
            # Edge widths based on score
            edge_weights = [subgraph[u][v].get("weight", 400) / 1000 * 2 for u, v in subgraph.edges()]
            
            # Draw network
            nodes = nx.draw_networkx_nodes(
                subgraph, pos, ax=ax,
                node_size=node_sizes,
                node_color=node_colors,
                cmap=plt.cm.YlOrRd,
                alpha=0.85
            )
            
            nx.draw_networkx_edges(
                subgraph, pos, ax=ax,
                width=edge_weights,
                alpha=0.4,
                edge_color='#95a5a6'
            )
            
            # Label top hub genes only
            top_labels = {n: n for n in list(subgraph.nodes())[:20] 
                         if n in (hub_genes.head(20)["gene"].tolist() if not hub_genes.empty else [])}
            
            if not top_labels:
                top_labels = {n: n for n in list(subgraph.nodes())[:15]}
            
            nx.draw_networkx_labels(
                subgraph, pos, labels=top_labels, ax=ax,
                font_size=8, font_weight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7)
            )
            
            plt.colorbar(nodes, ax=ax, label="Hub Score", fraction=0.046)
            ax.set_title(f"PPI Network - {disease_name}\n(Top {subgraph.number_of_nodes()} genes, {subgraph.number_of_edges()} interactions)",
                        fontsize=14, fontweight='bold')
            ax.axis('off')
            
            plt.tight_layout()
            plot_path = os.path.join(self.output_dir, f"ppi_network_{disease_name.replace(' ', '_')}.png")
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"  🕸️  PPI network saved: {plot_path}")
            
        except Exception as e:
            print(f"  ⚠️  Network visualization error: {e}")

    def analyze_network(self, gene_list: list, disease_name: str) -> pd.DataFrame:
        """Full PPI network analysis pipeline."""
        # Get interactions
        interactions = self.get_string_interactions(gene_list)
        
        if interactions.empty:
            print("  ⚠️  No PPI data available, returning genes as-is")
            return pd.DataFrame({"gene": gene_list[:20], "hub_score": [0.5] * min(20, len(gene_list))})
        
        # Build network
        G = self.build_network(interactions)
        
        if G.number_of_nodes() == 0:
            return pd.DataFrame()
        
        # Find hub genes
        hub_genes = self.calculate_hub_genes(G)
        
        # Visualize
        self.visualize_network(G, hub_genes, disease_name)
        
        return hub_genes
