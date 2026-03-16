"""
AI Clinical Biomarker Finder
Main Pipeline Orchestrator

Usage:
    python main.py --disease "diabetes mellitus" --datasets 3
    python main.py --disease "coronary artery disease"
    python main.py --disease "diabetic retinopathy" --datasets 2 --top-biomarkers 15
"""

import argparse
import os
import sys
import time
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Add modules to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from geo_retriever import GEODatasetRetriever
from deg_analyzer import DifferentialExpressionAnalyzer
from enrichment_analyzer import FunctionalEnrichmentAnalyzer
from ppi_analyzer import PPINetworkAnalyzer
from biomarker_scorer import BiomarkerScorer
from literature_miner import LiteratureMiner
from report_generator import ReportGenerator


def print_banner():
    banner = """
╔══════════════════════════════════════════════════════════════╗
║                                                              ║
║        🧬  AI CLINICAL BIOMARKER FINDER  🧬                 ║
║                                                              ║
║   GEO Database → DEG Analysis → PPI Network → DeepSeek R1   ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
"""
    print(banner)


def validate_config():
    """Validate required environment variables."""
    ncbi_email = os.getenv("NCBI_EMAIL")
    ollama_host = os.getenv("OLLAMA_HOST", "http://localhost:11434")
    ollama_model = os.getenv("OLLAMA_MODEL", "deepseek-r1")

    missing = []
    if not ncbi_email or ncbi_email == "your_email@example.com":
        missing.append("NCBI_EMAIL")

    if missing:
        print(f"❌ Missing configuration: {', '.join(missing)}")
        print("  Please edit your .env file with valid credentials.")
        print("  See .env.example for reference.")
        sys.exit(1)

    print(f"✅ Configuration loaded")
    print(f"   Email:        {ncbi_email}")
    print(f"   Ollama Host:  {ollama_host}")
    print(f"   Ollama Model: {ollama_model}")


def run_pipeline(disease_name: str,
                  n_datasets: int = 3,
                  top_biomarkers: int = 20,
                  pvalue_thresh: float = 0.05,
                  logfc_thresh: float = 1.0):
    """
    Run the complete biomarker discovery pipeline.

    Args:
        disease_name: Disease to analyze (e.g., "diabetes mellitus")
        n_datasets: Number of GEO datasets to download and analyze
        top_biomarkers: Number of top biomarkers to report
        pvalue_thresh: Adjusted p-value threshold for DEGs
        logfc_thresh: log2 fold change threshold for DEGs
    """

    print_banner()

    # Load configuration
    ncbi_email = os.getenv("NCBI_EMAIL")
    ncbi_api_key = os.getenv("NCBI_API_KEY")

    # Ollama configuration (replaces Anthropic)
    ollama_host = os.getenv("OLLAMA_HOST", "http://localhost:11434")
    ollama_model = os.getenv("OLLAMA_MODEL", "deepseek-r1")

    # Directories
    output_dir = os.getenv("OUTPUT_DIR", "./outputs")
    reports_dir = os.getenv("REPORTS_DIR", "./reports")
    data_dir = os.getenv("DATA_DIR", "./data")

    for d in [output_dir, reports_dir, data_dir]:
        os.makedirs(d, exist_ok=True)

    print(f"\n🎯 Target Disease: {disease_name}")
    print(f"📋 Settings: {n_datasets} datasets, p<{pvalue_thresh}, |logFC|>{logfc_thresh}")
    print("=" * 60)

    start_time = time.time()

    # ──────────────────────────────────────────────
    # STEP 1: Retrieve Gene Expression Datasets
    # ──────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("STEP 1/7: Dataset Retrieval from NCBI GEO")
    print("─" * 60)

    retriever = GEODatasetRetriever(email=ncbi_email, api_key=ncbi_api_key)
    datasets = retriever.retrieve_datasets_for_disease(disease_name, n_datasets=n_datasets)

    if not datasets:
        print("\n❌ No datasets retrieved. Possible causes:")
        print("   - Disease name too specific (try simpler name)")
        print("   - Network connectivity issue")
        print("   - NCBI servers unavailable")
        sys.exit(1)

    # ──────────────────────────────────────────────
    # STEP 2: Differential Gene Expression Analysis
    # ──────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("STEP 2/7: Differential Gene Expression Analysis")
    print("─" * 60)

    deg_analyzer = DifferentialExpressionAnalyzer(output_dir=output_dir)
    all_degs = []

    for dataset in datasets:
        degs = deg_analyzer.analyze_dataset(
            dataset,
            pvalue_threshold=pvalue_thresh,
            logfc_threshold=logfc_thresh
        )
        if not degs.empty:
            all_degs.append(degs)

    if not all_degs:
        print("⚠️  No significant DEGs found with current thresholds.")
        print("   Try relaxing thresholds: --pvalue 0.1 --logfc 0.5")
        # Use top genes by p-value even if not significant
        all_degs_relaxed = []
        for dataset in datasets:
            degs = deg_analyzer.analyze_dataset(dataset, pvalue_threshold=0.5, logfc_threshold=0.3)
            if not degs.empty:
                all_degs_relaxed.append(degs)
        all_degs = all_degs_relaxed

    # Find shared DEGs across datasets
    min_datasets_threshold = max(1, len(all_degs) - 1)
    shared_degs = deg_analyzer.find_shared_degs(all_degs, min_datasets=min_datasets_threshold)

    if shared_degs.empty and all_degs:
        shared_degs = all_degs[0]

    print(f"\n  📊 Total candidate genes: {len(shared_degs)}")

    # Generate heatmap
    if not shared_degs.empty:
        deg_analyzer.generate_heatmap(datasets, shared_degs, disease_name)

    # Get gene symbols for downstream analysis
    if "gene_id" in shared_degs.columns:
        gene_list = shared_degs["gene_id"].tolist()
    else:
        gene_list = shared_degs.iloc[:, 0].tolist()

    gene_list = [str(g) for g in gene_list if g][:200]  # Limit size

    # ──────────────────────────────────────────────
    # STEP 3: Functional Enrichment Analysis
    # ──────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("STEP 3/7: Functional Enrichment Analysis (GO, KEGG)")
    print("─" * 60)

    enrichment_analyzer = FunctionalEnrichmentAnalyzer(output_dir=output_dir)
    enrichment_results = enrichment_analyzer.run_enrichr_analysis(
        gene_list[:100],  # Enrichr works best with 100-200 genes
        description=f"{disease_name.replace(' ', '_')}_DEGs"
    )

    if enrichment_results:
        enrichment_analyzer.generate_enrichment_plots(enrichment_results, disease_name)

    top_pathways = enrichment_analyzer.get_top_pathways_summary(enrichment_results)

    # ──────────────────────────────────────────────
    # STEP 4: PPI Network Analysis
    # ──────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("STEP 4/7: Protein-Protein Interaction Network (STRING)")
    print("─" * 60)

    ppi_analyzer = PPINetworkAnalyzer(output_dir=output_dir)
    hub_genes = ppi_analyzer.analyze_network(gene_list, disease_name)

    # ──────────────────────────────────────────────
    # STEP 5: Biomarker Scoring
    # ──────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("STEP 5/7: Biomarker Scoring & Ranking")
    print("─" * 60)

    scorer = BiomarkerScorer(output_dir=output_dir)
    biomarkers_df = scorer.compute_biomarker_scores(
        degs=shared_degs,
        hub_genes=hub_genes,
        enrichment_results=enrichment_results,
        top_n=top_biomarkers
    )

    if biomarkers_df.empty:
        print("❌ Could not compute biomarker scores")
        sys.exit(1)

    # Print quick summary
    print("\n" + scorer.format_biomarker_panel(biomarkers_df, disease_name))

    # ──────────────────────────────────────────────
    # STEP 6: LLM Literature Analysis (DeepSeek R1 via Ollama)
    # ──────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("STEP 6/7: AI Literature Mining (DeepSeek R1 via Ollama + PubMed)")
    print("─" * 60)

    miner = LiteratureMiner(
        ollama_host=ollama_host,
        ollama_model=ollama_model,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key
    )

    # Convert biomarkers to list of dicts
    biomarker_list = biomarkers_df.to_dict(orient="records")

    # Analyze top 10 genes with LLM
    literature_analyses = miner.analyze_top_biomarkers(
        biomarker_list,
        disease_name,
        top_n=10
    )

    # Generate disease overview
    disease_overview = miner.generate_disease_overview(
        disease_name,
        biomarker_list,
        top_pathways
    )

    # ──────────────────────────────────────────────
    # STEP 7: Generate Reports
    # ──────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("STEP 7/7: Generating Reports")
    print("─" * 60)

    reporter = ReportGenerator(output_dir=output_dir, reports_dir=reports_dir)

    # HTML Report
    html_path = reporter.generate_html_report(
        disease_name=disease_name,
        datasets_info=datasets,
        biomarkers_df=biomarkers_df,
        enrichment_results=enrichment_results,
        literature_analyses=literature_analyses,
        disease_overview=disease_overview,
        hub_genes=hub_genes
    )

    # JSON Results
    json_path = reporter.save_results_json(
        disease_name=disease_name,
        biomarkers_df=biomarkers_df,
        enrichment_results=enrichment_results,
        literature_analyses=literature_analyses
    )

    # CSV Biomarkers
    csv_path = reporter.save_biomarkers_csv(biomarkers_df, disease_name)

    # ──────────────────────────────────────────────
    # Final Summary
    # ──────────────────────────────────────────────
    elapsed = time.time() - start_time

    print("\n" + "=" * 60)
    print("✅  ANALYSIS COMPLETE!")
    print("=" * 60)
    print(f"⏱️  Total time: {elapsed/60:.1f} minutes")
    print(f"\n📁 Output Files:")
    print(f"   🌐 HTML Report:  {html_path}")
    print(f"   📊 CSV Results:  {csv_path}")
    print(f"   📄 JSON Data:    {json_path}")

    # try to open report automatically in default browser
    try:
        import webbrowser
        webbrowser.open(f"file://{os.path.abspath(html_path)}")
    except Exception:
        pass

    # helpful message when user attempted to browse to 0.0.0.0:8000
    print("\n📌 Note: the pipeline does not start a web server by default.")
    print("   To view the HTML report you can either open it directly with your browser")
    print("   (e.g. file://<path>/" + os.path.basename(html_path) + ") or run:")
    print("      cd " + reports_dir + " && python -m http.server 8000")
    print("   then navigate to http://localhost:8000/" + os.path.basename(html_path))
    print(f"\n🧬 Top 10 Biomarker Candidates for {disease_name}:")

    for _, row in biomarkers_df.head(10).iterrows():
        arrow = "↑" if row.get("direction", "UP") == "UP" else "↓"
        print(f"   {arrow} {row['gene']:<12} Score: {row['biomarker_score']:.3f}")

    print("\n⚠️  Disclaimer: Results require experimental validation before clinical use.")
    print("=" * 60)

    return {
        "biomarkers": biomarkers_df,
        "html_report": html_path,
        "csv_path": csv_path,
        "json_path": json_path
    }


def main():
    parser = argparse.ArgumentParser(
        description="AI Clinical Biomarker Finder",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py --disease "diabetes mellitus"
  python main.py --disease "coronary artery disease" --datasets 3
  python main.py --disease "diabetic retinopathy" --top-biomarkers 15
  python main.py --disease "breast cancer" --pvalue 0.05 --logfc 1.5
        """
    )

    parser.add_argument(
        "--disease", "-d",
        type=str,
        required=True,
        help="Disease name to analyze (e.g., 'diabetes mellitus')"
    )
    parser.add_argument(
        "--datasets", "-n",
        type=int,
        default=3,
        help="Number of GEO datasets to download (default: 3)"
    )
    parser.add_argument(
        "--top-biomarkers", "-t",
        type=int,
        default=20,
        help="Number of top biomarkers to report (default: 20)"
    )
    parser.add_argument(
        "--pvalue", "-p",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold for DEGs (default: 0.05)"
    )
    parser.add_argument(
        "--logfc", "-l",
        type=float,
        default=1.0,
        help="log2 Fold Change threshold (default: 1.0)"
    )

    args = parser.parse_args()

    # Validate config
    validate_config()

    # Run pipeline
    run_pipeline(
        disease_name=args.disease,
        n_datasets=args.datasets,
        top_biomarkers=args.top_biomarkers,
        pvalue_thresh=args.pvalue,
        logfc_thresh=args.logfc
    )


if __name__ == "__main__":
    main()