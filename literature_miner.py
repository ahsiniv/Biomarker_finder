"""
Module 6: LLM Literature Evidence Mining
Uses Claude AI to search PubMed and analyze biomarker evidence.
"""

import requests
import json
import os
import time
import xml.etree.ElementTree as ET
from typing import Optional


class LiteratureMiner:
    def __init__(self, ollama_host: str, ollama_model: str, ncbi_email: str, ncbi_api_key: str = None):
        """
        Initialize literature miner.
        
        Args:
            ollama_host: Ollama server URL
            ollama_model: Model name (e.g., deepseek-r2)
            ncbi_email: Email for NCBI queries
            ncbi_api_key: Optional NCBI API key
        """
        self.ollama_host = ollama_host.rstrip('/')
        self.ollama_model = ollama_model
        self.ncbi_email = ncbi_email
        self.ncbi_api_key = ncbi_api_key
        self.pubmed_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def search_pubmed(self, query: str, max_results: int = 5) -> list:
        """
        Search PubMed for abstracts about a gene and disease.
        
        Args:
            query: Search query (e.g., "VEGFA diabetic retinopathy biomarker")
            max_results: Maximum results to return
            
        Returns:
            List of abstract dicts
        """
        # Search for PMIDs
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "sort": "relevance",
            "email": self.ncbi_email,
            "retmode": "json"
        }
        
        if self.ncbi_api_key:
            params["api_key"] = self.ncbi_api_key
        
        try:
            search_response = requests.get(
                f"{self.pubmed_base}/esearch.fcgi",
                params=params,
                timeout=15
            )
            search_data = search_response.json()
            pmids = search_data.get("esearchresult", {}).get("idlist", [])
            
            if not pmids:
                return []
            
            # Fetch abstracts
            fetch_params = {
                "db": "pubmed",
                "id": ",".join(pmids),
                "rettype": "abstract",
                "retmode": "xml",
                "email": self.ncbi_email
            }
            
            if self.ncbi_api_key:
                fetch_params["api_key"] = self.ncbi_api_key
            
            time.sleep(0.5)
            fetch_response = requests.get(
                f"{self.pubmed_base}/efetch.fcgi",
                params=fetch_params,
                timeout=20
            )
            
            # Parse XML abstracts
            abstracts = self._parse_pubmed_xml(fetch_response.text)
            return abstracts
            
        except Exception as e:
            return []

    def _parse_pubmed_xml(self, xml_text: str) -> list:
        """Parse PubMed XML response to extract abstracts."""
        articles = []
        
        try:
            root = ET.fromstring(xml_text)
            
            for article in root.findall(".//PubmedArticle"):
                try:
                    # Title
                    title_el = article.find(".//ArticleTitle")
                    title = title_el.text if title_el is not None else ""
                    
                    # Abstract
                    abstract_texts = article.findall(".//AbstractText")
                    abstract = " ".join([
                        el.text for el in abstract_texts if el.text
                    ])
                    
                    # Year
                    year_el = article.find(".//PubDate/Year")
                    year = year_el.text if year_el is not None else "N/A"
                    
                    # Journal
                    journal_el = article.find(".//Journal/Title")
                    journal = journal_el.text if journal_el is not None else ""
                    
                    if title and abstract:
                        articles.append({
                            "title": title,
                            "abstract": abstract[:1000],  # Limit length
                            "year": year,
                            "journal": journal
                        })
                        
                except Exception:
                    continue
                    
        except Exception:
            pass
        
        return articles

    def analyze_biomarker_with_llm(self, gene: str, disease: str, abstracts: list, deg_info: dict = None) -> dict:
        """
        Use Ollama to analyze biomarker evidence from literature.
        """
        # Prepare context
        abstract_text = ""
        for i, abs_data in enumerate(abstracts[:3], 1):
            abstract_text += f"\n[Study {i}] ({abs_data['year']}) {abs_data['title']}\n"
            abstract_text += f"Abstract: {abs_data['abstract']}\n"

        deg_text = ""
        if deg_info:
            deg_text = f"""
Computational Evidence:
- log2 Fold Change: {deg_info.get('logFC', 'N/A')}
- Adjusted p-value: {deg_info.get('adj_pvalue', 'N/A')}
- Direction: {deg_info.get('direction', 'N/A')} (in disease vs control)
- Biomarker Score: {deg_info.get('biomarker_score', 'N/A')}
"""

        prompt = f"""You are a clinical bioinformatics expert. Analyze the biomarker potential of {gene} in {disease}.

{deg_text}

Literature Evidence:
{abstract_text if abstract_text else "No specific literature found. Use your knowledge."}

Please provide a structured analysis with these exact sections:

1. GENE FUNCTION: Brief description of {gene}'s normal biological role (2-3 sentences)

2. DISEASE MECHANISM: How {gene} is involved in {disease} pathophysiology (2-3 sentences)

3. BIOMARKER POTENTIAL: Assessment of {gene} as a clinical biomarker for {disease} (2-3 sentences)

4. CLINICAL DETECTABILITY: Can this be measured clinically? (blood, urine, tissue?) (1-2 sentences)

5. EVIDENCE STRENGTH: Rate as HIGH/MEDIUM/LOW with brief justification

6. KEY PATHWAYS: List 2-3 pathways this gene is involved in

Be concise and evidence-based. If upregulated in disease, explain why high expression indicates pathology."""

        try:
            response = requests.post(
                f"{self.ollama_host}/api/generate",
                json={
                    "model": self.ollama_model,
                    "prompt": prompt,
                    "stream": False
                },
                timeout=120
            )
            response.raise_for_status()
            result = response.json()
            analysis_text = result.get("response", "Analysis unavailable.")

            # Parse evidence strength
            strength = "MEDIUM"
            if "EVIDENCE STRENGTH: HIGH" in analysis_text.upper():
                strength = "HIGH"
            elif "EVIDENCE STRENGTH: LOW" in analysis_text.upper():
                strength = "LOW"

            return {
                "gene": gene,
                "disease": disease,
                "analysis": analysis_text,
                "evidence_strength": strength,
                "n_pubmed_hits": len(abstracts)
            }
        except Exception as e:
            return {
                "gene": gene,
                "disease": disease,
                "analysis": f"Analysis unavailable: {e}",
                "evidence_strength": "UNKNOWN",
                "n_pubmed_hits": 0
            }

    def analyze_top_biomarkers(self, biomarkers: list, disease_name: str, top_n: int = 10) -> list:
        """
        Analyze top biomarkers using LLM + PubMed.
        
        Args:
            biomarkers: List of dicts with gene info
            disease_name: Disease name
            top_n: Number of biomarkers to analyze
            
        Returns:
            List of analysis dicts
        """
        print(f"\n🤖 AI Literature Analysis (Claude + PubMed)...")
        analyses = []
        
        for bm in biomarkers[:top_n]:
            gene = bm.get("gene", "")
            if not gene:
                continue
            
            print(f"  📚 Analyzing {gene}...")
            
            # Search PubMed
            query = f"{gene} {disease_name} biomarker"
            abstracts = self.search_pubmed(query, max_results=5)
            
            if not abstracts:
                query = f"{gene} {disease_name}"
                abstracts = self.search_pubmed(query, max_results=3)
            
            # Analyze with LLM
            analysis = self.analyze_biomarker_with_llm(gene, disease_name, abstracts, bm)
            analyses.append(analysis)
            
            time.sleep(1)  # Rate limiting
        
        return analyses

    def generate_disease_overview(self, disease_name: str, top_biomarkers: list, enrichment_summary: list) -> str:
        """
        Generate overall disease biomarker summary using Ollama.
        """
        print(f"\n🤖 Generating disease overview with Ollama...")

        # Prepare biomarker list
        bm_text = "\n".join([
            f"- {bm['gene']}: score={bm.get('biomarker_score', 0):.3f}, "
            f"logFC={bm.get('logFC', 0):.2f}, direction={bm.get('direction', 'N/A')}"
            for bm in top_biomarkers[:15]
        ])

        # Prepare pathway list
        path_text = "\n".join([
            f"- {p['term']} (p={p['adj_pvalue']:.4f})"
            for p in enrichment_summary[:10]
        ]) if enrichment_summary else "No pathway data available"

        prompt = f"""You are a clinical genomics expert. Based on computational analysis of gene expression data, \
provide a comprehensive biomarker discovery summary for {disease_name}.

Top Candidate Biomarkers (ranked by composite score):
{bm_text}

Top Enriched Pathways:
{path_text}

Please write a structured research summary with these sections:

## Disease Overview
Brief overview of {disease_name} and its clinical burden (3-4 sentences)

## Key Findings
Interpret the top biomarker candidates and what their expression changes suggest about {disease_name} pathophysiology (4-5 sentences)

## Top Biomarker Panel
Recommend the best 5-7 gene panel for clinical use, explaining why each was selected

## Pathway Insights
What do the enriched pathways tell us about the disease mechanisms? (3-4 sentences)

## Clinical Applications
How could these biomarkers be used clinically? (diagnosis, prognosis, drug targets) (3-4 sentences)

## Limitations & Next Steps
What validation would be needed to translate these to clinical use? (2-3 sentences)

Write in professional scientific style suitable for a research paper."""

        try:
            response = requests.post(
                f"{self.ollama_host}/api/generate",
                json={
                    "model": self.ollama_model,
                    "prompt": prompt,
                    "stream": False
                },
                timeout=180
            )
            response.raise_for_status()
            result = response.json()
            return result.get("response", "Overview generation failed.")
        except Exception as e:
            return f"Overview generation failed: {e}"
