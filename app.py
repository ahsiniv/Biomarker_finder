"""
FastAPI Backend for AI Clinical Biomarker Finder
Provides REST API + SSE streaming for the pipeline.
"""

import asyncio
import json
import os
import sys
import uuid
import threading
import queue
import time
from contextlib import redirect_stdout
from datetime import datetime
import logging

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from dotenv import load_dotenv

load_dotenv()

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

app = FastAPI(title="AI Clinical Biomarker Finder", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files (outputs/reports directories)
os.makedirs("./outputs", exist_ok=True)
os.makedirs("./reports", exist_ok=True)
os.makedirs("./data", exist_ok=True)

app.mount("/outputs", StaticFiles(directory="./outputs"), name="outputs")
app.mount("/reports", StaticFiles(directory="./reports"), name="reports")

# ─── In-memory job store ───────────────────────────────────────────────────────
jobs: dict[str, dict] = {}
job_queues: dict[str, queue.Queue] = {}


class JobRequest(BaseModel):
    disease: str
    datasets: int = 3
    top_biomarkers: int = 20
    pvalue: float = 0.05
    logfc: float = 1.0


class PipelineCapture:
    """Captures stdout and sends events to the SSE queue."""

    def __init__(self, job_id: str, q: queue.Queue):
        self.job_id = job_id
        self.q = q
        self._buffer = ""
        self._out = sys.__stdout__
        self.logger = logging.getLogger(__name__)

    def write(self, text: str):
        # Mirror output to the real stdout for visibility
        self._out.write(text)

        self._buffer += text
        if "\n" in self._buffer:
            lines = self._buffer.split("\n")
            for line in lines[:-1]:
                if line.strip():
                    self.logger.debug(line)
                    self.q.put({"type": "log", "message": line})
            self._buffer = lines[-1]

    def flush(self):
        self._out.flush()
        if self._buffer.strip():
            self.logger.debug(self._buffer)
            self.q.put({"type": "log", "message": self._buffer})
            self._buffer = ""

def run_pipeline_thread(job_id: str, req: JobRequest):
    """Run the pipeline in a background thread."""
    q = job_queues[job_id]

    def send(event_type: str, data: dict):
        q.put({"type": event_type, **data})

    jobs[job_id]["status"] = "running"
    jobs[job_id]["started_at"] = datetime.now().isoformat()

    capture = PipelineCapture(job_id, q)

    try:
        with redirect_stdout(capture):
            ncbi_email = os.getenv("NCBI_EMAIL")
            ncbi_api_key = os.getenv("NCBI_API_KEY")
            ollama_host = os.getenv("OLLAMA_HOST", "http://localhost:11434")
            ollama_model = os.getenv("OLLAMA_MODEL", "deepseek-r1")
            output_dir = os.getenv("OUTPUT_DIR", "./outputs")
            reports_dir = os.getenv("REPORTS_DIR", "./reports")

            send("step", {"step": 1, "label": "Retrieving GEO Datasets"})
            from geo_retriever import GEODatasetRetriever
            retriever = GEODatasetRetriever(email=ncbi_email, api_key=ncbi_api_key)
            datasets = retriever.retrieve_datasets_for_disease(req.disease, n_datasets=req.datasets)

            if not datasets:
                raise Exception("No GEO datasets found. Try a simpler disease name.")

            jobs[job_id]["datasets_count"] = len(datasets)
            send("progress", {"value": 15, "label": f"Found {len(datasets)} datasets"})

            send("step", {"step": 2, "label": "Differential Expression Analysis"})
            from deg_analyzer import DifferentialExpressionAnalyzer
            deg_analyzer = DifferentialExpressionAnalyzer(output_dir=output_dir)
            all_degs = []
            for dataset in datasets:
                degs = deg_analyzer.analyze_dataset(dataset, pvalue_threshold=req.pvalue, logfc_threshold=req.logfc)
                if not degs.empty:
                    all_degs.append(degs)

            if not all_degs:
                all_degs_relaxed = []
                for dataset in datasets:
                    degs = deg_analyzer.analyze_dataset(dataset, pvalue_threshold=0.5, logfc_threshold=0.3)
                    if not degs.empty:
                        all_degs_relaxed.append(degs)
                all_degs = all_degs_relaxed

            min_thresh = max(1, len(all_degs) - 1)
            shared_degs = deg_analyzer.find_shared_degs(all_degs, min_datasets=min_thresh)
            if shared_degs.empty and all_degs:
                shared_degs = all_degs[0]

            if not shared_degs.empty:
                deg_analyzer.generate_heatmap(datasets, shared_degs, req.disease)

            gene_list = shared_degs["gene_id"].tolist() if "gene_id" in shared_degs.columns else shared_degs.iloc[:, 0].tolist()
            gene_list = [str(g) for g in gene_list if g][:200]

            send("progress", {"value": 35, "label": f"{len(gene_list)} DEGs identified"})

            send("step", {"step": 3, "label": "Enrichment Analysis (GO/KEGG)"})
            from enrichment_analyzer import FunctionalEnrichmentAnalyzer
            enrichment_analyzer = FunctionalEnrichmentAnalyzer(output_dir=output_dir)
            enrichment_results = enrichment_analyzer.run_enrichr_analysis(
                gene_list[:100],
                description=f"{req.disease.replace(' ', '_')}_DEGs"
            )
            if enrichment_results:
                enrichment_analyzer.generate_enrichment_plots(enrichment_results, req.disease)
            top_pathways = enrichment_analyzer.get_top_pathways_summary(enrichment_results)

            send("progress", {"value": 55, "label": "Enrichment analysis complete"})

            send("step", {"step": 4, "label": "PPI Network Analysis (STRING)"})
            from ppi_analyzer import PPINetworkAnalyzer
            ppi_analyzer = PPINetworkAnalyzer(output_dir=output_dir)
            hub_genes = ppi_analyzer.analyze_network(gene_list, req.disease)

            send("progress", {"value": 70, "label": "Network built"})

            send("step", {"step": 5, "label": "Biomarker Scoring"})
            from biomarker_scorer import BiomarkerScorer
            scorer = BiomarkerScorer(output_dir=output_dir)
            biomarkers_df = scorer.compute_biomarker_scores(
                degs=shared_degs,
                hub_genes=hub_genes,
                enrichment_results=enrichment_results,
                top_n=req.top_biomarkers
            )

            if biomarkers_df.empty:
                raise Exception("Could not compute biomarker scores.")

            send("progress", {"value": 80, "label": f"Top {len(biomarkers_df)} biomarkers scored"})

            send("step", {"step": 6, "label": "AI Literature Mining"})
            from literature_miner import LiteratureMiner
            miner = LiteratureMiner(
                ollama_host=ollama_host, ollama_model=ollama_model,
                ncbi_email=ncbi_email, ncbi_api_key=ncbi_api_key
            )
            biomarker_list = biomarkers_df.to_dict(orient="records")
            literature_analyses = miner.analyze_top_biomarkers(biomarker_list, req.disease, top_n=10)
            disease_overview = miner.generate_disease_overview(req.disease, biomarker_list, top_pathways)

            send("progress", {"value": 90, "label": "Literature analysis complete"})

            send("step", {"step": 7, "label": "Generating Reports"})
            from report_generator import ReportGenerator
            reporter = ReportGenerator(output_dir=output_dir, reports_dir=reports_dir)
            html_path = reporter.generate_html_report(
                disease_name=req.disease, datasets_info=datasets, biomarkers_df=biomarkers_df,
                enrichment_results=enrichment_results, literature_analyses=literature_analyses,
                disease_overview=disease_overview, hub_genes=hub_genes
            )
            json_path = reporter.save_results_json(req.disease, biomarkers_df, enrichment_results, literature_analyses)
            csv_path = reporter.save_biomarkers_csv(biomarkers_df, req.disease)

            # Build results payload
            safe_disease = req.disease.replace(" ", "_")
            results = {
                "biomarkers": biomarker_list[:20],
                "enrichment": {cat: v[:10] for cat, v in enrichment_results.items()},
                "hub_genes": hub_genes.head(20).to_dict(orient="records") if not hub_genes.empty else [],
                "html_report": f"/reports/biomarker_report_{safe_disease}.html",
                "csv_path": f"/reports/biomarkers_{safe_disease}.csv",
                "json_path": f"/reports/results_{safe_disease}.json",
                # Image paths
                "plots": {
                    "biomarker_scores": "/outputs/biomarker_scores.png",
                    "ppi_network": f"/outputs/ppi_network_{safe_disease}.png",
                    "heatmap": f"/outputs/heatmap_{safe_disease}.png",
                }
            }

            jobs[job_id]["results"] = results
            jobs[job_id]["status"] = "completed"
            jobs[job_id]["completed_at"] = datetime.now().isoformat()

            send("progress", {"value": 100, "label": "Analysis complete!"})
            send("complete", {"results": results})

    except Exception as e:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = str(e)
        q.put({"type": "error", "message": str(e)})
    finally:
        q.put(None)  # Signal end of stream

        # Cleanup completed/failed jobs after a retention period to avoid memory growth
        def _cleanup_job():
            jobs.pop(job_id, None)
            job_queues.pop(job_id, None)

        cleanup_delay = int(os.getenv("JOB_RETENTION_SECONDS", "3600"))
        t = threading.Timer(cleanup_delay, _cleanup_job)
        t.daemon = True
        t.start()


@app.get("/", response_class=HTMLResponse)
async def index():
    with open("frontend.html", "r", encoding="utf-8") as f:
        return f.read()


@app.post("/api/jobs")
async def create_job(req: JobRequest):
    job_id = str(uuid.uuid4())[:8]
    jobs[job_id] = {
        "id": job_id,
        "disease": req.disease,
        "status": "queued",
        "created_at": datetime.now().isoformat(),
        "settings": req.dict()
    }
    job_queues[job_id] = queue.Queue()

    thread = threading.Thread(target=run_pipeline_thread, args=(job_id, req), daemon=True)
    thread.start()

    return {"job_id": job_id}


@app.get("/api/jobs/{job_id}/stream")
async def stream_job(job_id: str):
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    q = job_queues.get(job_id)

    async def event_generator():
        while True:
            try:
                # Non-blocking get with async sleep
                try:
                    item = q.get_nowait()
                except queue.Empty:
                    await asyncio.sleep(0.2)
                    continue

                if item is None:
                    yield f"data: {json.dumps({'type': 'done'})}\n\n"
                    break

                yield f"data: {json.dumps(item)}\n\n"

            except Exception as e:
                yield f"data: {json.dumps({'type': 'error', 'message': str(e)})}\n\n"
                break

    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "X-Accel-Buffering": "no",
        }
    )


@app.get("/api/jobs/{job_id}")
async def get_job(job_id: str):
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    return jobs[job_id]


@app.get("/api/jobs")
async def list_jobs():
    return list(jobs.values())


@app.delete("/api/jobs/{job_id}")
async def delete_job(job_id: str):
    jobs.pop(job_id, None)
    job_queues.pop(job_id, None)
    return {"ok": True}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app:app", host="0.0.0.0", port=8001, reload=False)