# 🧬 AI Clinical Biomarker Finder

An end-to-end pipeline that takes a **disease name as input** and outputs **ranked clinical biomarker candidates** by integrating NCBI GEO gene expression data, differential expression analysis, pathway enrichment, PPI network analysis, and an LLM-based literature mining step (Ollama).

---

## 📁 Project Structure

```
biomarker_finder/
├── main.py                          # ← Main entry point (run this)
├── requirements.txt                 # Python dependencies
├── env.example                      # Config template
├── .env                             # Your config (create this)
│
├── modules/
│   ├── geo_retriever.py             # Step 1: Fetch GEO datasets
│   ├── deg_analyzer.py              # Step 2: Differential expression
│   ├── enrichment_analyzer.py       # Step 3: GO/KEGG enrichment (Enrichr API)
│   ├── ppi_analyzer.py              # Step 4: STRING PPI network
│   ├── biomarker_scorer.py          # Step 5: Composite scoring
│   ├── literature_miner.py          # Step 6: LLM literature mining (Ollama) + PubMed
│   └── report_generator.py          # Step 7: HTML/CSV/JSON reports
│
├── outputs/                         # Plots (volcano, heatmap, network)
├── reports/                         # Final reports (HTML, CSV, JSON)
└── data/                            # GEO dataset cache
```

---

## ⚙️ Setup Instructions

### 1. Install Python (3.10+)
Download from: https://python.org

### 2. Open in VS Code
```
File → Open Folder → Select biomarker_finder/
```

### 3. Create Virtual Environment
Open terminal in VS Code (`Ctrl+`` `):
```bash
python -m venv .venv
.venv\Scripts\activate          # Windows
source .venv/bin/activate       # Mac/Linux
```

### 4. Install Dependencies
```bash
pip install -r requirements.txt
```

### 5. Configure API Keys
```bash
copy env.example .env          # Windows
cp env.example .env            # Mac/Linux
```

Then edit `.env` with your credentials:

```
# Required for NCBI Entrez API usage
NCBI_EMAIL=you@email.com        # Your email (required by NCBI)
NCBI_API_KEY=...                # Optional: get from ncbi.nlm.nih.gov/account

# Ollama LLM (local or remote)
OLLAMA_HOST=http://localhost:11434
OLLAMA_MODEL=deepseek-r1
```

---

## 🚀 Running the Pipeline

### Via Terminal (recommended)
```bash
# Basic usage
python main.py --disease "diabetes mellitus"

# With options
python main.py --disease "coronary artery disease" --datasets 3

# More options
python main.py --disease "diabetic retinopathy" --datasets 2 --top-biomarkers 15 --pvalue 0.05 --logfc 1.0
```

### Via VS Code Debug
1. Press `F5` or go to Run → Start Debugging
2. Select a launch configuration (pre-configured for 3 diseases)
3. Watch the terminal output

### Command Line Arguments
| Argument | Default | Description |
|----------|---------|-------------|
| `--disease` | required | Disease name (e.g., "diabetes mellitus") |
| `--datasets` | 3 | Number of GEO datasets to download |
| `--top-biomarkers` | 20 | Number of biomarkers in output |
| `--pvalue` | 0.05 | Adjusted p-value threshold |
| `--logfc` | 1.0 | log2 Fold Change threshold |

---

## 📊 Output Files

After running, you'll find:

| File | Description |
|------|-------------|
| `reports/biomarker_report_[disease].html` | **Main report** - open in browser |
| `reports/biomarkers_[disease].csv` | Ranked biomarker list |
| `reports/results_[disease].json` | All data in JSON format |
| `outputs/volcano_*.png` | Volcano plots per dataset |
| `outputs/heatmap_*.png` | DEG heatmap |
| `outputs/ppi_network_*.png` | STRING network visualization |
| `outputs/enrichment_*.png` | Pathway enrichment plots |

---

## 🔬 Pipeline Overview

```
User Input: Disease Name
        ↓
Step 1: Search NCBI GEO → Download 3 datasets (GSE files)
        ↓
Step 2: Differential Expression Analysis
        → Welch's t-test + Benjamini-Hochberg FDR correction
        → Find genes with |logFC| > 1 AND adj.p < 0.05
        → Intersection across multiple datasets
        ↓
Step 3: Enrichment Analysis (Enrichr API)
        → GO Biological Process, Molecular Function
        → KEGG Pathways, Reactome, DisGeNET
        ↓
Step 4: PPI Network (STRING API)
        → Build protein interaction network
        → Calculate: Degree, Betweenness, PageRank, MCC
        → Identify hub genes
        ↓
Step 5: Biomarker Scoring
        → DEG significance (40%) + Network score (35%) + Pathway count (25%)
        ↓
Step 6: AI Literature Mining
        → Search PubMed for each top gene
        → Claude AI generates mechanistic analysis
        → Evidence strength assessment
        ↓
Step 7: Generate HTML Report + CSV + JSON
```

---

## 🧪 Example Diseases to Try
- `"diabetes mellitus"`
- `"diabetic retinopathy"`
- `"coronary artery disease"`
- `"breast cancer"`
- `"Alzheimer disease"`
- `"rheumatoid arthritis"`
- `"hypertension"`
- `"liver fibrosis"`

---

## 📚 Data Sources Used
- **NCBI GEO** - Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/)
- **STRING** - Protein interaction database (https://string-db.org)
- **Enrichr** - Enrichment analysis (https://maayanlab.cloud/Enrichr)
- **PubMed** - Literature search (https://pubmed.ncbi.nlm.nih.gov)
- **Ollama LLM** - Literature synthesis (https://ollama.ai)

---

## ⚠️ Important Notes
- Analysis time: ~15-30 minutes depending on disease and network speed
- GEO downloads are cached in `data/geo_cache/` - rerunning the same disease is faster
- Results are for **research purposes only** - require experimental validation
- NCBI requires a valid email address for API access

---

## 🔧 Troubleshooting

**"No datasets found"**
→ Try broader disease name: `"diabetes"` instead of `"type 2 diabetes mellitus"`

**"ModuleNotFoundError"**
→ Make sure venv is activated and `pip install -r requirements.txt` was run

**"Invalid API key"**  
→ Check your `.env` file has the correct `OLLAMA_HOST` and `OLLAMA_MODEL` values (or ensure your Ollama server is running).

**Analysis takes too long**
→ Reduce `--datasets 1` for faster results
