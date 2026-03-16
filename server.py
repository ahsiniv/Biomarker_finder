import asyncio
import json
import os
from uuid import uuid4
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import StreamingResponse, FileResponse, JSONResponse
from typing import Dict, Any

from main import run_pipeline

app = FastAPI()

# job store: job_id -> dict
jobs: Dict[str, Dict[str, Any]] = {}

# helpers for event streaming

def make_event_queue():
    return asyncio.Queue()


def enqueue_event(job_id: str, event: Dict[str, Any]):
    q = jobs[job_id].get("events")
    if q:
        q.put_nowait(event)


async def run_pipeline_job(job_id: str, payload: Dict[str, Any]):
    """Run the pipeline in a background thread, streaming logs back."""
    enqueue_event(job_id, {"type": "log", "message": "Starting job"})
    jobs[job_id]["status"] = "running"

    # create a writer that intercepts print() output
    class InterceptWriter:
        def write(self, s):
            # called with chunks; split into lines
            for line in s.splitlines():
                if not line:
                    continue
                enqueue_event(job_id, {"type": "log", "message": line})
                # detect step
                import re
                m = re.search(r"STEP\s*(\d+)/7", line)
                if m:
                    step_num = int(m.group(1))
                    label = line.split(":", 1)[-1].strip() if ":" in line else line
                    enqueue_event(job_id, {"type": "step", "step": step_num, "label": label})
                    prog = int((step_num / 7) * 100)
                    enqueue_event(job_id, {"type": "progress", "value": prog, "label": label})
        def flush(self):
            pass

    # run pipeline in thread, capturing stdout via redirect
    loop = asyncio.get_event_loop()
    try:
        def target():
            from contextlib import redirect_stdout
            buf = InterceptWriter()
            with redirect_stdout(buf):
                return run_pipeline(
                    disease_name=payload.get("disease"),
                    n_datasets=payload.get("datasets", 3),
                    top_biomarkers=payload.get("top_biomarkers", 20),
                    pvalue_thresh=payload.get("pvalue", 0.05),
                    logfc_thresh=payload.get("logfc", 1.0)
                )
        results = await loop.run_in_executor(None, target)
        # convert pandas objects to serializable forms
        try:
            import pandas as pd
            if isinstance(results.get("biomarkers"), pd.DataFrame):
                results["biomarkers"] = results["biomarkers"].to_dict(orient="records")
        except ImportError:
            pass
        jobs[job_id]["results"] = results
        enqueue_event(job_id, {"type": "complete", "results": results})
    except Exception as e:
        enqueue_event(job_id, {"type": "error", "message": str(e)})
        jobs[job_id]["status"] = "failed"
        return
    jobs[job_id]["status"] = "completed"
    enqueue_event(job_id, {"type": "done"})


@app.post("/api/jobs")
async def create_job(request: Request):
    payload = await request.json()
    job_id = str(uuid4())
    jobs[job_id] = {
        "status": "queued",
        "disease": payload.get("disease"),
        "results": None,
        "events": make_event_queue()
    }
    # start background task
    asyncio.create_task(run_pipeline_job(job_id, payload))
    return {"job_id": job_id}


@app.get("/api/jobs/{job_id}")
async def get_job(job_id: str):
    job = jobs.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return {"status": job["status"], "disease": job.get("disease"), "results": job.get("results")}


@app.get("/api/jobs/{job_id}/stream")
async def stream(job_id: str):
    job = jobs.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    q = job["events"]

    async def event_generator():
        try:
            while True:
                event = await q.get()
                yield f"data: {json.dumps(event)}\n\n"
                if event.get("type") in ("complete", "error", "done"):
                    break
        except asyncio.CancelledError:
            pass

    return StreamingResponse(event_generator(), media_type="text/event-stream")


@app.get("/")
async def root():
    # serve the static frontend
    return FileResponse("frontend.html")


@app.get("/static/{path:path}")
async def static_files(path: str):
    # serve any other static file if needed
    return FileResponse(path)
