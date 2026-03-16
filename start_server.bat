@echo off
echo.
echo ╔══════════════════════════════════════════════════╗
echo ║   🧬  AI Clinical Biomarker Finder  - Web UI     ║
echo ╚══════════════════════════════════════════════════╝
echo.

:: Activate venv if present
if exist ".venv\Scripts\activate.bat" (
    call .venv\Scripts\activate.bat
    echo Virtual environment activated
)

:: Ensure .env exists
if not exist ".env" (
    if exist "env.example" (
        copy "env.example" ".env" >nul
        echo ⚠️  Created .env from env.example — please update it with your credentials
    ) else (
        echo ❌ No .env file found. Create one from env.example.
        exit /b 1
    )
)

:: Install dependencies
echo Installing dependencies (requirements.txt)...
pip install -q -r requirements.txt

echo.
echo Starting server at http://localhost:8000
echo Press Ctrl+C to stop
echo.

uvicorn app:app --host 0.0.0.0 --port 8000 --reload
pause
