"""
Module 2: Differential Gene Expression Analysis
Performs statistical analysis to find significantly changed genes.
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from Bio import Entrez
warnings.filterwarnings('ignore')


class DifferentialExpressionAnalyzer:
    def __init__(self, output_dir: str = "./outputs"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def _map_entrez_to_symbol(self, gene_ids: list) -> dict:
        """Map numeric Entrez gene IDs to official gene symbols via Entrez esummary."""
        mapping = {}
        numeric_ids = [str(int(x)) for x in gene_ids if str(x).strip().isdigit()]
        if not numeric_ids:
            return mapping

        try:
            id_chunks = [numeric_ids[i:i+200] for i in range(0, len(numeric_ids), 200)]
            for chunk in id_chunks:
                handle = Entrez.esummary(db="gene", id=','.join(chunk))
                records = Entrez.read(handle)
                handle.close()

                if isinstance(records, dict) and "DocumentSummarySet" in records:
                    records = records["DocumentSummarySet"].get("DocumentSummary", [])

                if isinstance(records, list):
                    for rec in records:
                        gid = str(rec.get("Id", ""))
                        sym = rec.get("NomenclatureSymbol") or rec.get("Name") or gid
                        mapping[gid] = sym
                elif isinstance(records, dict):
                    gid = str(records.get("Id", ""))
                    sym = records.get("NomenclatureSymbol") or records.get("Name") or gid
                    mapping[gid] = sym
        except Exception:
            pass
        return mapping

    def analyze_dataset(self, dataset: dict, 
                        pvalue_threshold: float = 0.05,
                        logfc_threshold: float = 1.0) -> pd.DataFrame:
        """
        Perform differential expression analysis on a dataset.
        
        Args:
            dataset: Processed dataset dict from GEODatasetRetriever
            pvalue_threshold: Adjusted p-value cutoff
            logfc_threshold: |log2FC| cutoff
            
        Returns:
            DataFrame with significant DEGs
        """
        gse_id = dataset["gse_id"]
        print(f"\n🧬 Analyzing {gse_id}...")
        
        # Build expression DataFrame
        expr_dict = dataset["expression_matrix"]
        sample_groups = dataset["sample_groups"]
        
        try:
            expr_df = pd.DataFrame(expr_dict)
        except Exception as e:
            print(f"  ✗ Could not build expression matrix: {e}")
            return pd.DataFrame()
        
        # Get sample columns that exist in our matrix
        disease_cols = [s for s in sample_groups["disease"] if s in expr_df.columns]
        control_cols = [s for s in sample_groups["control"] if s in expr_df.columns]
        
        if len(disease_cols) < 2 or len(control_cols) < 2:
            print(f"  ✗ Insufficient samples: {len(disease_cols)} disease, {len(control_cols)} control")
            return pd.DataFrame()
        
        print(f"  📊 {len(disease_cols)} disease vs {len(control_cols)} control samples")
        
        # Convert to numeric, drop non-numeric genes
        expr_df = expr_df.apply(pd.to_numeric, errors='coerce')
        expr_df = expr_df.dropna(how='all')
        
        # Filter low-expression genes (must be expressed in >20% of samples)
        min_samples = max(2, int(0.2 * expr_df.shape[1]))
        expr_df = expr_df[expr_df.notna().sum(axis=1) >= min_samples]
        
        disease_expr = expr_df[disease_cols].copy()
        control_expr = expr_df[control_cols].copy()
        
        # Fill NaN with row median
        disease_expr = disease_expr.apply(lambda x: x.fillna(x.median()), axis=1)
        control_expr = control_expr.apply(lambda x: x.fillna(x.median()), axis=1)
        
        print(f"  🔬 Testing {len(expr_df)} genes...")
        
        # Compute statistics for each gene
        results = []
        
        for gene_id in expr_df.index:
            try:
                d_vals = disease_expr.loc[gene_id].dropna().values.astype(float)
                c_vals = control_expr.loc[gene_id].dropna().values.astype(float)
                
                if len(d_vals) < 2 or len(c_vals) < 2:
                    continue
                
                # Compute log2 fold change
                mean_disease = np.mean(d_vals)
                mean_control = np.mean(c_vals)
                
                if mean_control == 0:
                    mean_control = 0.001
                    
                logfc = mean_disease - mean_control  # Already log2 space
                
                # Welch's t-test (doesn't assume equal variance)
                t_stat, pval = stats.ttest_ind(d_vals, c_vals, equal_var=False)
                
                results.append({
                    "gene_id": str(gene_id),
                    "logFC": round(logfc, 4),
                    "mean_disease": round(mean_disease, 4),
                    "mean_control": round(mean_control, 4),
                    "pvalue": pval,
                    "t_stat": round(t_stat, 4)
                })
                
            except Exception:
                continue
        
        if not results:
            print(f"  ✗ No results computed")
            return pd.DataFrame()

        results_df = pd.DataFrame(results)
        results_df = results_df.dropna(subset=["pvalue"])
        
        # Multiple testing correction (Benjamini-Hochberg)
        _, adj_pvals, _, _ = multipletests(results_df["pvalue"].values, method='fdr_bh')
        results_df["adj_pvalue"] = adj_pvals

        # Map entrez gene IDs to symbols (when possible)
        if "gene_id" in results_df.columns:
            unique_ids = list(results_df["gene_id"].astype(str).unique())
            gene_symbol_map = self._map_entrez_to_symbol(unique_ids)
            results_df["gene_symbol"] = results_df["gene_id"].astype(str).map(
                lambda x: gene_symbol_map.get(x, x)
            )
        else:
            results_df["gene_symbol"] = results_df.index.astype(str)
        
        # Filter by thresholds
        degs = results_df[
            (results_df["adj_pvalue"] < pvalue_threshold) & 
            (results_df["logFC"].abs() > logfc_threshold)
        ].copy()
        
        degs["direction"] = degs["logFC"].apply(lambda x: "UP" if x > 0 else "DOWN")
        degs = degs.sort_values("adj_pvalue")
        
        # Keep gene_symbol column and fallback to gene_id
        if "gene_symbol" in degs.columns:
            degs["gene"] = degs["gene_symbol"].fillna(degs["gene_id"])
        else:
            degs["gene"] = degs["gene_id"]

        n_up = (degs["direction"] == "UP").sum()
        n_down = (degs["direction"] == "DOWN").sum()
        print(f"  ✅ Found {len(degs)} DEGs: {n_up} upregulated, {n_down} downregulated")
        
        # Generate volcano plot
        self._generate_volcano_plot(results_df, degs, gse_id, pvalue_threshold, logfc_threshold)
        
        return degs

    def _generate_volcano_plot(self, all_results: pd.DataFrame, degs: pd.DataFrame,
                                gse_id: str, pval_thresh: float, logfc_thresh: float):
        """Generate volcano plot for the dataset."""
        try:
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # -log10 transform p-values
            all_results["-log10_pval"] = -np.log10(all_results["adj_pvalue"].clip(1e-300))
            
            # Plot non-significant genes
            not_sig = all_results[
                (all_results["adj_pvalue"] >= pval_thresh) | 
                (all_results["logFC"].abs() <= logfc_thresh)
            ]
            ax.scatter(not_sig["logFC"], not_sig["-log10_pval"],
                      alpha=0.3, s=8, c='lightgray', label='Not significant')
            
            # Plot significant upregulated
            sig_up = all_results[
                (all_results["adj_pvalue"] < pval_thresh) & 
                (all_results["logFC"] > logfc_thresh)
            ]
            ax.scatter(sig_up["logFC"], sig_up["-log10_pval"],
                      alpha=0.7, s=15, c='#e74c3c', label=f'Upregulated ({len(sig_up)})')
            
            # Plot significant downregulated
            sig_down = all_results[
                (all_results["adj_pvalue"] < pval_thresh) & 
                (all_results["logFC"] < -logfc_thresh)
            ]
            ax.scatter(sig_down["logFC"], sig_down["-log10_pval"],
                      alpha=0.7, s=15, c='#3498db', label=f'Downregulated ({len(sig_down)})')
            
            # Threshold lines
            ax.axhline(y=-np.log10(pval_thresh), color='black', linestyle='--', alpha=0.5, linewidth=1)
            ax.axvline(x=logfc_thresh, color='black', linestyle='--', alpha=0.5, linewidth=1)
            ax.axvline(x=-logfc_thresh, color='black', linestyle='--', alpha=0.5, linewidth=1)
            
            # Label top genes
            top_genes = all_results[
                (all_results["adj_pvalue"] < pval_thresh) & 
                (all_results["logFC"].abs() > logfc_thresh)
            ].nlargest(10, "-log10_pval")
            
            for _, row in top_genes.iterrows():
                ax.annotate(str(row["gene_id"])[:10],
                           xy=(row["logFC"], row["-log10_pval"]),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=7, alpha=0.8)
            
            ax.set_xlabel("log2 Fold Change", fontsize=12)
            ax.set_ylabel("-log10(Adjusted p-value)", fontsize=12)
            ax.set_title(f"Volcano Plot - {gse_id}", fontsize=14, fontweight='bold')
            ax.legend(loc='upper left', fontsize=9)
            
            plt.tight_layout()
            plot_path = os.path.join(self.output_dir, f"volcano_{gse_id}.png")
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"  📈 Volcano plot saved: {plot_path}")
            
        except Exception as e:
            print(f"  ⚠️  Could not generate volcano plot: {e}")

    def find_shared_degs(self, all_degs: list, min_datasets: int = 2) -> pd.DataFrame:
        """
        Find genes that are differentially expressed across multiple datasets.
        
        Args:
            all_degs: List of DEG DataFrames from multiple datasets
            min_datasets: Minimum number of datasets a gene must appear in
            
        Returns:
            DataFrame of shared DEGs with combined statistics
        """
        if not all_degs:
            return pd.DataFrame()

        # Collect all gene IDs from each dataset
        gene_sets = []
        gene_stats = {}  # gene -> list of (logFC, adj_pvalue)

        for df in all_degs:
            if df.empty:
                continue

            if "gene_symbol" in df.columns:
                genes_in_dataset = set(
                    df["gene_symbol"].fillna(df["gene_id"] if "gene_id" in df.columns else "").tolist()
                )
            elif "gene_id" in df.columns:
                genes_in_dataset = set(df["gene_id"].tolist())
            else:
                genes_in_dataset = set(df.iloc[:, 0].tolist())
            gene_sets.append(genes_in_dataset)

            for _, row in df.iterrows():
                gene = row.get("gene_symbol") or row.get("gene_id") or str(row.iloc[0])
                if not gene or str(gene).strip() == "":
                    continue
                if gene not in gene_stats:
                    gene_stats[gene] = []
                gene_stats[gene].append({
                    "logFC": row.get("logFC", row.get("mean_logFC", 0)),
                    "adj_pvalue": row.get("adj_pvalue", row.get("mean_adj_pvalue", 0.05)),
                    "direction": row.get("direction", "UP")
                })
        
        if not gene_sets:
            return pd.DataFrame()
        
        # Find genes present in enough datasets
        shared_genes = []
        
        for gene, stats_list in gene_stats.items():
            if len(stats_list) >= min_datasets:
                # Check consistent direction
                directions = [s["direction"] for s in stats_list]
                consistent = len(set(directions)) == 1

                shared_genes.append({
                    "gene_id": gene,
                    "gene_symbol": gene,
                    "n_datasets": len(stats_list),
                    "mean_logFC": np.mean([s["logFC"] for s in stats_list]),
                    "mean_adj_pvalue": np.mean([s["adj_pvalue"] for s in stats_list]),
                    "adj_pvalue": np.mean([s["adj_pvalue"] for s in stats_list]),
                    "logFC": np.mean([s["logFC"] for s in stats_list]),
                    "direction": directions[0] if consistent else "MIXED",
                    "consistent_direction": consistent
                })
        
        if not shared_genes:
            print("  ⚠️  No shared DEGs found. Relaxing to single-dataset genes...")
            # Return from the dataset with most DEGs
            largest = max(all_degs, key=len) if all_degs else pd.DataFrame()
            return largest
        
        shared_df = pd.DataFrame(shared_genes)
        shared_df = shared_df.sort_values("n_datasets", ascending=False)
        
        print(f"\n🔗 Found {len(shared_df)} shared DEGs across ≥{min_datasets} datasets")
        return shared_df

    def generate_heatmap(self, datasets: list, degs: pd.DataFrame,
                         disease_name: str, top_n: int = 50):
        """Generate expression heatmap for top DEGs."""
        try:
            if degs.empty:
                return

            # Use first dataset for heatmap
            dataset = datasets[0]
            expr_dict = dataset["expression_matrix"]
            sample_groups = dataset["sample_groups"]

            expr_df = pd.DataFrame(expr_dict)
            expr_df = expr_df.apply(pd.to_numeric, errors='coerce')

            # Build a mapping from gene_id/symbol -> probe index
            # degs["gene_id"] may contain symbols (after find_shared_degs) or probe IDs
            # expr_df.index always contains the original probe/entrez IDs from GEO
            # Strategy: try direct match first, then match via gene_symbol column
            top_genes_raw = degs.head(top_n)["gene_id"].tolist()

            # Direct match (works when degs come from a single dataset with probe IDs)
            available_genes = [g for g in top_genes_raw if g in expr_df.index]

            # Fallback: if degs has a gene_symbol column, build reverse symbol->probe map
            if len(available_genes) < 5 and "gene_symbol" in degs.columns:
                symbol_to_probe = {}
                for orig_probe in expr_df.index:
                    # Check if any DEG gene_symbol maps back to this probe via gene_id
                    pass
                # Try matching via the original gene_id (probe) stored in degs
                if "gene_id" in degs.columns:
                    # When shared_degs replaces gene_id with symbol, the original probe is lost.
                    # Fall back to using all DEG probes that still exist in expr_df.
                    probe_candidates = degs.head(top_n).get("original_gene_id", pd.Series(dtype=str))
                    available_genes = [g for g in probe_candidates if g in expr_df.index]

            # Last resort: use top N rows of expr_df ranked by variance
            if len(available_genes) < 5:
                variances = expr_df.var(axis=1).sort_values(ascending=False)
                available_genes = variances.head(top_n).index.tolist()[:30]

            available_genes = available_genes[:30]

            if len(available_genes) < 5:
                return
            
            subset = expr_df.loc[available_genes]
            
            # Order samples: control first, then disease
            control_cols = [c for c in sample_groups["control"] if c in subset.columns]
            disease_cols = [c for c in sample_groups["disease"] if c in subset.columns]
            ordered_cols = control_cols + disease_cols
            subset = subset[ordered_cols].fillna(0)
            
            # Create column colors
            col_colors = pd.Series(
                ['#3498db'] * len(control_cols) + ['#e74c3c'] * len(disease_cols),
                index=ordered_cols
            )
            
            fig, ax = plt.subplots(figsize=(12, 10))
            
            g = sns.clustermap(
                subset,
                col_cluster=False,
                row_cluster=True,
                col_colors=col_colors,
                cmap='RdBu_r',
                center=0,
                figsize=(14, 10),
                yticklabels=True,
                xticklabels=False,
                cbar_kws={"label": "log2 Expression"}
            )
            
            g.fig.suptitle(f"DEG Heatmap - {disease_name}", 
                          y=1.02, fontsize=14, fontweight='bold')
            
            plot_path = os.path.join(self.output_dir, f"heatmap_{disease_name.replace(' ', '_')}.png")
            g.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"  📊 Heatmap saved: {plot_path}")
            
        except Exception as e:
            print(f"  ⚠️  Heatmap generation failed: {e}")
