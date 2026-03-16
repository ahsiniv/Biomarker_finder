"""
Module 5: Biomarker Scoring System
Ranks candidate biomarkers using multiple evidence sources.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')


class BiomarkerScorer:
    def __init__(self, output_dir: str = "./outputs"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def compute_biomarker_scores(self, 
                                  degs: pd.DataFrame,
                                  hub_genes: pd.DataFrame,
                                  enrichment_results: dict,
                                  top_n: int = 20) -> pd.DataFrame:
        """
        Compute composite biomarker score for each candidate gene.
        
        Score components:
        - DEG significance (adjusted p-value, logFC)
        - Network centrality (hub score from STRING)
        - Pathway involvement (number of enriched pathways)
        
        Args:
            degs: Differentially expressed genes DataFrame
            hub_genes: Hub genes from PPI network
            enrichment_results: Enrichment analysis results
            top_n: Number of top biomarkers to return
            
        Returns:
            Ranked biomarker DataFrame
        """
        print(f"\n🏆 Computing biomarker scores...")
        
        if degs.empty:
            print("  ✗ No DEG data available")
            return pd.DataFrame()
        
        # --- Component 1: DEG Score ---
        deg_score_col = "gene_id" if "gene_id" in degs.columns else degs.columns[0]
        
        deg_scores = {}
        for _, row in degs.iterrows():
            gene = row.get("gene_id", row.iloc[0])
            
            # Score based on significance and fold change
            adj_pval = row.get("adj_pvalue", row.get("mean_adj_pvalue", 0.05))
            logfc = abs(row.get("mean_logFC", row.get("logFC", 1)))
            
            # -log10(p-value) normalized
            sig_score = min(-np.log10(max(adj_pval, 1e-50)), 50) / 50
            fc_score = min(logfc / 5, 1.0)  # Cap at logFC=5
            
            # Number of datasets gene appears in
            n_datasets = row.get("n_datasets", 1)
            dataset_score = min(n_datasets / 3, 1.0)
            
            deg_scores[gene] = {
                "sig_score": sig_score,
                "fc_score": fc_score,
                "dataset_score": dataset_score,
                "logFC": logfc,
                "adj_pvalue": adj_pval,
                "direction": row.get("direction", "UP"),
                "n_datasets": n_datasets
            }
        
        # --- Component 2: Network Score ---
        network_scores = {}
        if not hub_genes.empty and "gene" in hub_genes.columns:
            for _, row in hub_genes.iterrows():
                gene = row["gene"]
                network_scores[gene] = row.get("hub_score", 0)
        
        # --- Component 3: Pathway Score ---
        pathway_gene_counts = {}
        for category, results in enrichment_results.items():
            for pathway in results:
                genes_in_pathway = pathway.get("genes", [])
                for g in genes_in_pathway:
                    pathway_gene_counts[g] = pathway_gene_counts.get(g, 0) + 1
        
        max_pathway_count = max(pathway_gene_counts.values()) if pathway_gene_counts else 1
        
        # --- Combine Scores ---
        biomarker_list = []
        
        all_genes = set(deg_scores.keys())
        
        for gene in all_genes:
            deg_info = deg_scores.get(gene, {})
            
            # DEG component (40% weight)
            deg_component = (
                deg_info.get("sig_score", 0) * 0.5 +
                deg_info.get("fc_score", 0) * 0.3 +
                deg_info.get("dataset_score", 0) * 0.2
            ) * 0.40
            
            # Network component (35% weight)
            net_score = network_scores.get(gene, 0)
            network_component = net_score * 0.35
            
            # Pathway component (25% weight)
            pathway_count = pathway_gene_counts.get(gene, 0)
            pathway_component = (pathway_count / max_pathway_count) * 0.25
            
            # Total composite score
            total_score = deg_component + network_component + pathway_component
            
            biomarker_list.append({
                "gene": gene,
                "biomarker_score": round(total_score, 4),
                "deg_score": round(deg_component / 0.40, 4),
                "network_score": round(net_score, 4),
                "pathway_count": pathway_count,
                "logFC": deg_info.get("logFC", 0),
                "adj_pvalue": deg_info.get("adj_pvalue", 1.0),
                "direction": deg_info.get("direction", "UP"),
                "n_datasets": deg_info.get("n_datasets", 1)
            })
        
        if not biomarker_list:
            return pd.DataFrame()
        
        result_df = pd.DataFrame(biomarker_list)
        result_df = result_df.sort_values("biomarker_score", ascending=False)
        result_df = result_df.head(top_n)
        
        # Add rank
        result_df["rank"] = range(1, len(result_df) + 1)
        
        print(f"  ✅ Top 5 biomarker candidates:")
        for _, row in result_df.head(5).iterrows():
            print(f"     #{int(row['rank'])}. {row['gene']} (score: {row['biomarker_score']:.3f}, "
                  f"logFC: {row['logFC']:.2f}, dir: {row['direction']})")
        
        # Generate visualization
        self._plot_biomarker_scores(result_df)
        
        return result_df

    def _plot_biomarker_scores(self, df: pd.DataFrame):
        """Generate biomarker ranking visualization."""
        try:
            top = df.head(15)
            
            fig, axes = plt.subplots(1, 2, figsize=(16, 8))
            
            # Left: Composite score bar chart
            ax1 = axes[0]
            colors = ['#e74c3c' if d == 'UP' else '#3498db' for d in top['direction']]
            bars = ax1.barh(
                range(len(top)),
                top['biomarker_score'],
                color=colors,
                alpha=0.8,
                edgecolor='black',
                linewidth=0.5
            )
            ax1.set_yticks(range(len(top)))
            ax1.set_yticklabels(top['gene'], fontsize=10, fontweight='bold')
            ax1.set_xlabel("Composite Biomarker Score", fontsize=11)
            ax1.set_title("Biomarker Ranking\n(Red=Upregulated, Blue=Downregulated)",
                         fontsize=12, fontweight='bold')
            ax1.invert_yaxis()
            
            # Add score labels
            for bar, score in zip(bars, top['biomarker_score']):
                ax1.text(bar.get_width() + 0.005, bar.get_y() + bar.get_height()/2,
                        f'{score:.3f}', va='center', fontsize=8)
            
            # Right: Score components
            ax2 = axes[1]
            genes = top.head(10)['gene'].tolist()
            deg_s = top.head(10)['deg_score'].tolist()
            net_s = top.head(10)['network_score'].tolist()
            path_s = (top.head(10)['pathway_count'] / max(top['pathway_count'].max(), 1)).tolist()
            
            x = np.arange(len(genes))
            width = 0.25
            
            ax2.bar(x - width, deg_s, width, label='DEG Score', color='#e74c3c', alpha=0.8)
            ax2.bar(x, net_s, width, label='Network Score', color='#2ecc71', alpha=0.8)
            ax2.bar(x + width, path_s, width, label='Pathway Score', color='#9b59b6', alpha=0.8)
            
            ax2.set_xticks(x)
            ax2.set_xticklabels(genes, rotation=45, ha='right', fontsize=9)
            ax2.set_ylabel("Score Component", fontsize=11)
            ax2.set_title("Score Components Breakdown", fontsize=12, fontweight='bold')
            ax2.legend(fontsize=9)
            ax2.set_ylim(0, 1.1)
            
            plt.suptitle("Biomarker Analysis Summary", fontsize=14, fontweight='bold', y=1.02)
            plt.tight_layout()
            
            plot_path = os.path.join(self.output_dir, "biomarker_scores.png")
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"  📊 Biomarker score plot saved: {plot_path}")
            
        except Exception as e:
            print(f"  ⚠️  Score plot error: {e}")

    def format_biomarker_panel(self, df: pd.DataFrame, disease_name: str) -> str:
        """Format biomarker panel as readable text."""
        if df.empty:
            return "No biomarkers identified."
        
        lines = [
            f"PROPOSED BIOMARKER PANEL FOR {disease_name.upper()}",
            "=" * 60,
            ""
        ]
        
        for _, row in df.head(10).iterrows():
            direction_symbol = "↑" if row["direction"] == "UP" else "↓"
            lines.append(
                f"#{int(row['rank']):2d}. {row['gene']:<12} {direction_symbol}  "
                f"Score: {row['biomarker_score']:.3f}  |  "
                f"logFC: {row['logFC']:+.2f}  |  "
                f"adj.p: {row['adj_pvalue']:.2e}  |  "
                f"Pathways: {row['pathway_count']}"
            )
        
        lines.append("")
        lines.append("↑ = Upregulated in disease   ↓ = Downregulated in disease")
        
        return "\n".join(lines)
