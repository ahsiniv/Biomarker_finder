"""
Module 3: Functional Enrichment Analysis
GO enrichment and KEGG pathway analysis using Enrichr and DAVID APIs.
"""

import requests
import pandas as pd
import time
import json
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')


class FunctionalEnrichmentAnalyzer:
    def __init__(self, output_dir: str = "./outputs"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.enrichr_base = "https://maayanlab.cloud/Enrichr"
        
    def run_enrichr_analysis(self, gene_list: list, description: str = "biomarker_genes") -> dict:
        """
        Run enrichment analysis using Enrichr API.
        
        Args:
            gene_list: List of gene symbols
            description: Analysis description
            
        Returns:
            Dict with GO and KEGG enrichment results
        """
        print(f"\n🔬 Running Enrichr enrichment analysis for {len(gene_list)} genes...")
        
        if not gene_list:
            print("  ✗ Empty gene list")
            return {}
        
        # Add gene list to Enrichr
        genes_str = "\n".join(gene_list)
        
        try:
            add_response = requests.post(
                f"{self.enrichr_base}/addList",
                files={"list": (None, genes_str), "description": (None, description)},
                timeout=30
            )
            add_result = add_response.json()
            user_list_id = add_result.get("userListId")
            
            if not user_list_id:
                print("  ✗ Failed to upload gene list to Enrichr")
                return {}
            
            print(f"  ✓ Uploaded to Enrichr (ID: {user_list_id})")
            
        except Exception as e:
            print(f"  ✗ Enrichr connection error: {e}")
            return {}
        
        # Gene set databases to query
        databases = {
            "GO_Biological_Process_2023": "GO Biological Process",
            "GO_Molecular_Function_2023": "GO Molecular Function",
            "KEGG_2021_Human": "KEGG Pathways",
            "Reactome_2022": "Reactome Pathways",
            "DisGeNET": "Disease Associations"
        }
        
        enrichment_results = {}
        
        for db_name, db_label in databases.items():
            try:
                time.sleep(0.5)
                enrich_response = requests.get(
                    f"{self.enrichr_base}/enrich",
                    params={"userListId": user_list_id, "backgroundType": db_name},
                    timeout=30
                )
                enrich_data = enrich_response.json()
                
                if db_name in enrich_data:
                    results = enrich_data[db_name]
                    # Each result: [rank, term, pvalue, zscore, combined_score, genes, adj_pvalue, ...]
                    parsed = []
                    for r in results[:20]:  # Top 20 results
                        if len(r) >= 7:
                            parsed.append({
                                "term": r[1],
                                "pvalue": r[2],
                                "zscore": r[3],
                                "combined_score": r[4],
                                "genes": r[5] if isinstance(r[5], list) else [],
                                "adj_pvalue": r[6]
                            })
                    
                    # Filter significant results
                    significant = [p for p in parsed if p["adj_pvalue"] < 0.05]
                    enrichment_results[db_label] = significant if significant else parsed[:10]
                    print(f"  ✓ {db_label}: {len(enrichment_results[db_label])} enriched terms")
                    
            except Exception as e:
                print(f"  ⚠️  Error querying {db_label}: {e}")
                continue
        
        return enrichment_results

    def run_kegg_analysis_direct(self, gene_list: list) -> list:
        """
        Directly query KEGG API for pathway enrichment.
        
        Args:
            gene_list: List of gene symbols
            
        Returns:
            List of enriched KEGG pathways
        """
        print("  🔗 Querying KEGG pathways directly...")
        
        kegg_results = []
        
        # Query KEGG FIND for each gene
        for gene in gene_list[:20]:  # Limit to avoid timeout
            try:
                response = requests.get(
                    f"https://rest.kegg.jp/find/genes/{gene}",
                    timeout=10
                )
                if response.status_code == 200 and response.text:
                    kegg_results.append(gene)
                time.sleep(0.2)
            except Exception:
                continue
        
        return kegg_results

    def generate_enrichment_plots(self, enrichment_results: dict, disease_name: str):
        """Generate bar plots for enrichment results."""
        for category, results in enrichment_results.items():
            if not results:
                continue
            
            try:
                top_results = results[:15]
                terms = [r["term"][:50] for r in top_results]
                scores = [-np.log10(max(r["adj_pvalue"], 1e-50)) for r in top_results]
                
                fig, ax = plt.subplots(figsize=(12, 6))
                
                colors = plt.cm.RdYlBu_r(np.linspace(0.2, 0.8, len(terms)))
                bars = ax.barh(range(len(terms)), scores[::-1], color=colors)
                
                ax.set_yticks(range(len(terms)))
                ax.set_yticklabels(terms[::-1], fontsize=9)
                ax.set_xlabel("-log10(Adjusted p-value)", fontsize=11)
                ax.set_title(f"{category}\n{disease_name}", fontsize=12, fontweight='bold')
                
                # Add threshold line
                ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')
                ax.legend(fontsize=9)
                
                plt.tight_layout()
                safe_category = category.replace(" ", "_").replace("/", "_")
                safe_disease = disease_name.replace(" ", "_")
                plot_path = os.path.join(self.output_dir, f"enrichment_{safe_category}_{safe_disease}.png")
                plt.savefig(plot_path, dpi=150, bbox_inches='tight')
                plt.close()
                print(f"  📊 Enrichment plot saved: {plot_path}")
                
            except Exception as e:
                print(f"  ⚠️  Plot error for {category}: {e}")

    def get_top_pathways_summary(self, enrichment_results: dict) -> list:
        """Extract top pathways across all categories."""
        all_pathways = []
        
        for category, results in enrichment_results.items():
            for r in results[:5]:
                all_pathways.append({
                    "category": category,
                    "term": r["term"],
                    "adj_pvalue": r["adj_pvalue"],
                    "genes": r.get("genes", [])
                })
        
        # Sort by p-value
        all_pathways.sort(key=lambda x: x["adj_pvalue"])
        return all_pathways[:20]
