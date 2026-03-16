#!/bin/bash
# ─────────────────────────────────────────────────────
#  AI Clinical Biomarker Finder — Web Server Startup
# ─────────────────────────────────────────────────────

echo ""
echo "╔══════════════════════════════════════════════════╗"
echo "║   🧬  AI Clinical Biomarker Finder  — Web UI     ║"
echo "╚══════════════════════════════════════════════════╝"
echo ""

# Check for .env
if [ ! -f ".env" ]; then
  if [ -f "env.example" ]; then
    cp env.example .env
    echo "⚠️  Created .env from env.example — please update it with your credentials"
  else
    echo "❌ No .env file found. Create one from env.example."
    exit 1
  fi
fi

# Activate venv if present
if [ -d ".venv" ]; then
  source .venv/bin/activate 2>/dev/null || source .venv/Scripts/activate 2>/dev/null
  echo "✅ Virtual environment activated"
fi

# Install / update deps
echo "📦 Installing dependencies (requirements.txt)..."
pip install -q -r requirements.txt

echo ""
echo "🌐 Starting server at http://localhost:8000"
echo "   Press Ctrl+C to stop"
echo ""

uvicorn app:app --host 0.0.0.0 --port 8000 --reload
