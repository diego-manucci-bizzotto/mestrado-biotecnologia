# Motif Scan

Sistema inicial para busca de motifs em promotores, com backend FastAPI e frontend React TS/Vite.

## Estrutura

```text
backend/   FastAPI, algoritmos de motif, p-value e q-value
frontend/  React TS/Vite com shadcn/ui
docs/      plano cientifico e referencias
```

## Backend

```powershell
cd backend
.\venv\Scripts\Activate.ps1
pip install -r requirements.txt
uvicorn app.main:app --reload
```

API: `http://127.0.0.1:8000/docs`

Se a porta `8000` ja estiver ocupada, use outra porta e aponte o frontend para ela:

```powershell
uvicorn app.main:app --reload --port 8010
```

## Frontend

```powershell
cd frontend
npm install
npm run dev
```

App: `http://localhost:5173`

Com backend em porta alternativa:

```powershell
$env:VITE_API_BASE_URL="http://127.0.0.1:8010/api"
npm run dev -- --port 5174
```
