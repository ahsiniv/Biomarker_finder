"""
Module 7: Report Generator
Generates comprehensive PDF and HTML reports.
"""

import os
import json
from datetime import datetime
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import warnings
warnings.filterwarnings('ignore')


class ReportGenerator:
    def __init__(self, output_dir: str = "./outputs", reports_dir: str = "./reports"):
        self.output_dir = output_dir
        self.reports_dir = reports_dir
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(reports_dir, exist_ok=True)

    def generate_html_report(self, 
                              disease_name: str,
                              datasets_info: list,
                              biomarkers_df: pd.DataFrame,
                              enrichment_results: dict,
                              literature_analyses: list,
                              disease_overview: str,
                              hub_genes: pd.DataFrame = None) -> str:
        """Generate comprehensive HTML report."""
        
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        safe_disease = disease_name.replace(" ", "_")
        
        # Build biomarker table HTML
        bm_rows = ""
        for _, row in biomarkers_df.head(20).iterrows():
            direction_badge = (
                f'<span style="background:#e74c3c;color:white;padding:2px 6px;border-radius:3px">↑ UP</span>'
                if row.get("direction", "UP") == "UP"
                else f'<span style="background:#3498db;color:white;padding:2px 6px;border-radius:3px">↓ DOWN</span>'
            )
            bm_rows += f"""
            <tr>
                <td><strong>#{int(row.get('rank', 0))}</strong></td>
                <td><strong style="color:#2c3e50">{row['gene']}</strong></td>
                <td>{row.get('biomarker_score', 0):.4f}</td>
                <td>{direction_badge}</td>
                <td>{row.get('logFC', 0):+.3f}</td>
                <td>{row.get('adj_pvalue', 1):.2e}</td>
                <td>{row.get('n_datasets', 1)}</td>
                <td>{int(row.get('pathway_count', 0))}</td>
                <td>{row.get('network_score', 0):.3f}</td>
            </tr>"""
        
        # Build enrichment table HTML
        enrichment_rows = ""
        for category, results in enrichment_results.items():
            for r in results[:5]:
                enrichment_rows += f"""
                <tr>
                    <td>{category}</td>
                    <td>{r['term'][:60] + ('...' if len(r['term']) > 60 else '')}</td>
                    <td>{r['adj_pvalue']:.2e}</td>
                    <td>{', '.join(r.get('genes', [])[:5])}</td>
                </tr>"""
        
        # Build literature section
        literature_html = ""
        for analysis in literature_analyses[:10]:
            strength_color = {
                "HIGH": "#27ae60",
                "MEDIUM": "#f39c12", 
                "LOW": "#e74c3c",
                "UNKNOWN": "#95a5a6"
            }.get(analysis.get("evidence_strength", "UNKNOWN"), "#95a5a6")
            
            lit_text = analysis.get("analysis", "").replace("\n", "<br>")
            
            literature_html += f"""
            <div class="gene-card">
                <div class="gene-header">
                    <h3>{analysis['gene']}</h3>
                    <span class="evidence-badge" style="background:{strength_color}">
                        {analysis.get('evidence_strength', 'N/A')} Evidence
                    </span>
                    <span class="pubmed-badge">📚 {analysis.get('n_pubmed_hits', 0)} PubMed studies</span>
                </div>
                <div class="gene-analysis">{lit_text}</div>
            </div>"""
        
        # Build datasets info
        datasets_html = ""
        for d in datasets_info:
            groups = d.get("sample_groups", {})
            datasets_html += f"""
            <div class="dataset-card">
                <strong>{d.get('gse_id', 'N/A')}</strong>
                <p>{d.get('title', '')[:80]}...</p>
                <p>Disease samples: {len(groups.get('disease', []))} | 
                   Control samples: {len(groups.get('control', []))}</p>
            </div>"""
        
        # Overview formatted
        overview_html = disease_overview.replace("\n## ", "<h3>").replace("\n\n", "</p><p>").replace("##", "<h3>").replace("\n", "<br>")
        
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AI Biomarker Report - {disease_name}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: 'Segoe UI', Arial, sans-serif; background: #f0f4f8; color: #2c3e50; }}
        
        .header {{ background: linear-gradient(135deg, #1a252f, #2c3e50, #34495e); 
                   color: white; padding: 50px 40px; text-align: center; }}
        .header h1 {{ font-size: 2.8em; margin-bottom: 10px; letter-spacing: -1px; }}
        .header h2 {{ font-size: 1.6em; color: #3498db; margin-bottom: 15px; }}
        .header .meta {{ font-size: 0.9em; color: #bdc3c7; }}
        .badge {{ display: inline-block; background: #3498db; color: white; 
                  padding: 5px 15px; border-radius: 20px; margin: 5px; font-size: 0.85em; }}
        
        .container {{ max-width: 1200px; margin: 0 auto; padding: 30px 20px; }}
        
        .section {{ background: white; border-radius: 12px; padding: 30px; margin-bottom: 25px;
                    box-shadow: 0 4px 15px rgba(0,0,0,0.08); }}
        .section h2 {{ font-size: 1.6em; color: #2c3e50; padding-bottom: 12px;
                       border-bottom: 3px solid #3498db; margin-bottom: 20px; }}
        
        table {{ width: 100%; border-collapse: collapse; font-size: 0.9em; }}
        th {{ background: #2c3e50; color: white; padding: 12px; text-align: left; }}
        td {{ padding: 10px 12px; border-bottom: 1px solid #ecf0f1; }}
        tr:hover {{ background: #f8f9fa; }}
        tr:nth-child(even) {{ background: #fdfdfe; }}
        
        .gene-card {{ border: 1px solid #e8ecef; border-radius: 10px; padding: 20px;
                      margin-bottom: 15px; transition: box-shadow 0.2s; }}
        .gene-card:hover {{ box-shadow: 0 4px 12px rgba(0,0,0,0.1); }}
        .gene-header {{ display: flex; align-items: center; gap: 12px; margin-bottom: 12px; }}
        .gene-header h3 {{ font-size: 1.3em; color: #2980b9; margin: 0; }}
        .evidence-badge {{ color: white; padding: 4px 10px; border-radius: 12px; font-size: 0.8em; font-weight: bold; }}
        .pubmed-badge {{ background: #ecf0f1; padding: 4px 10px; border-radius: 12px; font-size: 0.8em; }}
        .gene-analysis {{ font-size: 0.9em; line-height: 1.7; color: #555; }}
        
        .dataset-card {{ background: #f8f9fa; border-left: 4px solid #3498db;
                         padding: 12px 15px; margin-bottom: 10px; border-radius: 0 8px 8px 0; }}
        .dataset-card strong {{ color: #2980b9; font-size: 1.1em; }}
        .dataset-card p {{ margin-top: 5px; font-size: 0.9em; color: #555; }}
        
        .stats-grid {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 15px; margin-bottom: 20px; }}
        .stat-box {{ background: linear-gradient(135deg, #3498db, #2980b9); color: white;
                     border-radius: 10px; padding: 20px; text-align: center; }}
        .stat-box .number {{ font-size: 2.2em; font-weight: bold; }}
        .stat-box .label {{ font-size: 0.85em; opacity: 0.9; margin-top: 5px; }}
        
        .overview-content {{ line-height: 1.8; font-size: 0.95em; color: #444; }}
        .overview-content h3 {{ color: #2c3e50; margin: 20px 0 8px 0; font-size: 1.2em;
                                border-left: 4px solid #3498db; padding-left: 10px; }}
        
        .top-panel {{ background: linear-gradient(135deg, #f8f9fa, #ecf0f1);
                      border-radius: 10px; padding: 20px; margin-bottom: 20px; }}
        .top-panel h3 {{ color: #2c3e50; margin-bottom: 15px; }}
        .panel-gene {{ display: inline-flex; align-items: center; gap: 6px;
                       background: white; border: 2px solid #3498db; border-radius: 8px;
                       padding: 8px 14px; margin: 5px; font-weight: bold; }}
        
        .footer {{ background: #2c3e50; color: #bdc3c7; text-align: center; padding: 20px; margin-top: 30px; }}
    </style>
</head>
<body>

<div class="header">
    <h1>🧬 AI Clinical Biomarker Finder</h1>
    <h2>{disease_name}</h2>
    <div class="meta">
        Generated: {timestamp} | 
        Datasets analyzed: {len(datasets_info)} | 
        <span class="badge">GEO Database</span>
        <span class="badge">STRING PPI</span>
        <span class="badge">Enrichr</span>
        <span class="badge">Claude AI</span>
    </div>
</div>

<div class="container">

    <!-- Stats Summary -->
    <div class="section">
        <h2>📊 Analysis Summary</h2>
        <div class="stats-grid">
            <div class="stat-box">
                <div class="number">{len(datasets_info)}</div>
                <div class="label">GEO Datasets</div>
            </div>
            <div class="stat-box" style="background:linear-gradient(135deg,#e74c3c,#c0392b)">
                <div class="number">{len(biomarkers_df)}</div>
                <div class="label">Candidate Biomarkers</div>
            </div>
            <div class="stat-box" style="background:linear-gradient(135deg,#27ae60,#229954)">
                <div class="number">{sum(len(v) for v in enrichment_results.values())}</div>
                <div class="label">Enriched Pathways</div>
            </div>
            <div class="stat-box" style="background:linear-gradient(135deg,#9b59b6,#8e44ad)">
                <div class="number">{len(literature_analyses)}</div>
                <div class="label">Genes Literature-Reviewed</div>
            </div>
        </div>
    </div>

    <!-- Top Panel -->
    <div class="section">
        <h2>🏆 Proposed Biomarker Panel</h2>
        <div class="top-panel">
            <h3>Top Candidate Genes for {disease_name}</h3>
            {"".join([f'<span class="panel-gene">#{int(row["rank"])} {row["gene"]} {"↑" if row["direction"]=="UP" else "↓"}</span>' for _, row in biomarkers_df.head(8).iterrows()])}
        </div>
    </div>

    <!-- Datasets -->
    <div class="section">
        <h2>🗄️ Gene Expression Datasets</h2>
        {datasets_html}
    </div>

    <!-- Biomarker Rankings -->
    <div class="section">
        <h2>🧬 Biomarker Rankings</h2>
        <table>
            <thead>
                <tr>
                    <th>Rank</th><th>Gene</th><th>Score</th><th>Direction</th>
                    <th>log2FC</th><th>adj.p-value</th><th>Datasets</th>
                    <th>Pathways</th><th>Network</th>
                </tr>
            </thead>
            <tbody>
                {bm_rows}
            </tbody>
        </table>
    </div>

    <!-- Pathway Enrichment -->
    <div class="section">
        <h2>🔬 Pathway Enrichment Analysis</h2>
        <table>
            <thead>
                <tr><th>Category</th><th>Term</th><th>adj.p-value</th><th>Key Genes</th></tr>
            </thead>
            <tbody>
                {enrichment_rows}
            </tbody>
        </table>
    </div>

    <!-- Disease Overview -->
    <div class="section">
        <h2>🤖 AI-Generated Analysis Overview</h2>
        <div class="overview-content">
            {overview_html}
        </div>
    </div>

    <!-- Gene-level Literature Analysis -->
    <div class="section">
        <h2>📚 Gene-Level Literature Evidence</h2>
        {literature_html}
    </div>

</div>

<div class="footer">
    AI Clinical Biomarker Finder | Powered by NCBI GEO, STRING, Enrichr, Claude AI<br>
    <small>For research purposes only. Clinical validation required before medical use.</small>
</div>

</body>
</html>"""
        
        report_path = os.path.join(self.reports_dir, f"biomarker_report_{safe_disease}.html")
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(html)
        
        print(f"\n✅ HTML Report saved: {report_path}")
        return report_path

    def save_results_json(self, disease_name: str, biomarkers_df: pd.DataFrame,
                          enrichment_results: dict, literature_analyses: list) -> str:
        """Save all results as structured JSON."""
        safe_disease = disease_name.replace(" ", "_")
        
        results = {
            "disease": disease_name,
            "timestamp": datetime.now().isoformat(),
            "biomarker_panel": biomarkers_df.head(20).to_dict(orient="records") if not biomarkers_df.empty else [],
            "enrichment": {
                cat: results[:10] for cat, results in enrichment_results.items()
            },
            "literature": literature_analyses
        }
        
        json_path = os.path.join(self.reports_dir, f"results_{safe_disease}.json")
        with open(json_path, "w") as f:
            json.dump(results, f, indent=2, default=str)
        
        print(f"📄 JSON results saved: {json_path}")
        return json_path

    def save_biomarkers_csv(self, biomarkers_df: pd.DataFrame, disease_name: str) -> str:
        """Save biomarker rankings as CSV."""
        safe_disease = disease_name.replace(" ", "_")
        csv_path = os.path.join(self.reports_dir, f"biomarkers_{safe_disease}.csv")
        biomarkers_df.to_csv(csv_path, index=False)
        print(f"📊 CSV results saved: {csv_path}")
        return csv_path
