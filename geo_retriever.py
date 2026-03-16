"""
Module 1: GEO Dataset Retrieval
Fetches gene expression datasets from NCBI GEO based on disease name.
"""

import requests
import time
import os
import json
import GEOparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from Bio import Entrez
import xml.etree.ElementTree as ET


class GEODatasetRetriever:
    def __init__(self, email: str, api_key: str = None):
        """
        Initialize GEO retriever with NCBI credentials.
        
        Args:
            email: Your email (required by NCBI)
            api_key: NCBI API key (optional, increases rate limits)
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.cache_dir = "./data/geo_cache"
        os.makedirs(self.cache_dir, exist_ok=True)

    def search_geo_datasets(self, disease_name: str, max_results: int = 10) -> list:
        """
        Search NCBI GEO for datasets related to a disease.
        
        Args:
            disease_name: Name of the disease (e.g., 'diabetes mellitus')
            max_results: Maximum number of datasets to retrieve
            
        Returns:
            List of GSE accession IDs
        """
        print(f"\n🔍 Searching GEO for: '{disease_name}'...")
        
        # Search GEO datasets
        search_term = f'"{disease_name}"[Title/Abstract] AND "Homo sapiens"[Organism] AND "expression profiling by array"[DataSet Type]'
        
        try:
            handle = Entrez.esearch(db="gds", term=search_term, retmax=max_results * 2, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            geo_ids = record.get("IdList", [])
            
            if not geo_ids:
                # Fallback: broader search
                search_term = f'"{disease_name}"[Title/Abstract] AND "Homo sapiens"[Organism]'
                handle = Entrez.esearch(db="gds", term=search_term, retmax=max_results * 2, sort="relevance")
                record = Entrez.read(handle)
                handle.close()
                geo_ids = record.get("IdList", [])
            
            print(f"  ✓ Found {len(geo_ids)} datasets in GEO")
            return geo_ids[:max_results * 2]
            
        except Exception as e:
            print(f"  ✗ GEO search error: {e}")
            return []

    def fetch_dataset_metadata(self, geo_ids: list) -> list:
        """
        Fetch metadata for GEO datasets to filter relevant ones.
        
        Args:
            geo_ids: List of GEO numeric IDs
            
        Returns:
            List of dicts with GSE accessions and metadata
        """
        datasets = []
        
        if not geo_ids:
            return datasets
            
        print(f"\n📋 Fetching metadata for {len(geo_ids)} datasets...")
        
        for geo_id in tqdm(geo_ids[:15]):  # Limit to avoid rate limiting
            try:
                handle = Entrez.esummary(db="gds", id=geo_id)
                record = Entrez.read(handle)
                handle.close()
                
                if record:
                    summary = record[0]
                    accession = summary.get("Accession", "")
                    
                    # Only include GSE datasets (not GDS/GPL)
                    if accession.startswith("GSE"):
                        n_samples = int(summary.get("n_samples", 0))
                        
                        # Filter: must have at least 6 samples for meaningful analysis
                        if n_samples >= 6:
                            datasets.append({
                                "gse_id": accession,
                                "title": summary.get("title", ""),
                                "n_samples": n_samples,
                                "organism": summary.get("taxon", ""),
                                "summary": summary.get("summary", "")[:200]
                            })
                
                time.sleep(0.4)  # Be nice to NCBI servers
                
            except Exception as e:
                continue
        
        # Sort by sample size (more samples = more power)
        datasets.sort(key=lambda x: x["n_samples"], reverse=True)
        print(f"  ✓ {len(datasets)} valid datasets found with ≥6 samples")
        return datasets

    def download_and_parse_gse(self, gse_id: str) -> dict:
        """
        Download and parse a GSE dataset.
        
        Args:
            gse_id: GSE accession ID (e.g., 'GSE42148')
            
        Returns:
            Dict with expression matrix and sample groups
        """
        cache_file = os.path.join(self.cache_dir, f"{gse_id}_processed.json")
        
        # Check cache first
        if os.path.exists(cache_file):
            print(f"  📦 Loading cached data for {gse_id}...")
            with open(cache_file, 'r') as f:
                return json.load(f)
        
        print(f"  ⬇️  Downloading {gse_id}...")
        
        try:
            gse = GEOparse.get_GEO(
                geo=gse_id,
                destdir=self.cache_dir,
                silent=True
            )
            
            # Extract expression matrix
            expr_data = self._extract_expression_matrix(gse)
            
            if expr_data is None:
                return None
            
            # Identify control vs disease groups
            sample_groups = self._identify_sample_groups(gse)
            
            result = {
                "gse_id": gse_id,
                "title": gse.metadata.get("title", [""])[0],
                "expression_matrix": expr_data,
                "sample_groups": sample_groups,
                "n_genes": len(expr_data),
                "n_samples": sum(len(v) for v in sample_groups.values())
            }
            
            # Cache result
            with open(cache_file, 'w') as f:
                json.dump(result, f)
            
            return result
            
        except Exception as e:
            print(f"  ✗ Error downloading {gse_id}: {e}")
            return None

    def _extract_expression_matrix(self, gse) -> dict:
        """Extract and normalize expression matrix from GSE object."""
        try:
            pivot_samples = gse.pivot_samples("VALUE")
            
            if pivot_samples is None or pivot_samples.empty:
                return None
            
            # Remove rows with all NaN
            pivot_samples = pivot_samples.dropna(how='all')
            
            # Log2 transform if values look like they're not log-transformed
            if pivot_samples.max().max() > 100:
                pivot_samples = np.log2(pivot_samples + 1)
            
            return pivot_samples.to_dict()
            
        except Exception as e:
            return None

    def _identify_sample_groups(self, gse) -> dict:
        """Identify disease vs control sample groups from metadata."""
        disease_keywords = [
            'disease', 'patient', 'case', 'affected', 'diabetic', 
            'tumor', 'cancer', 'infected', 'treated', 'hypertension',
            'coronary', 'retinopathy', 'failure', 'disorder', 'syndrome'
        ]
        control_keywords = [
            'control', 'healthy', 'normal', 'wild', 'untreated', 
            'baseline', 'mock', 'vehicle', 'placebo'
        ]
        
        disease_samples = []
        control_samples = []
        
        for gsm_id, gsm in gse.gsms.items():
            title = gsm.metadata.get("title", [""])[0].lower()
            characteristics = " ".join([
                " ".join(v) for v in gsm.metadata.get("characteristics_ch1", [])
            ]).lower()
            source = " ".join(gsm.metadata.get("source_name_ch1", [""])).lower()
            
            full_text = f"{title} {characteristics} {source}"
            
            is_control = any(kw in full_text for kw in control_keywords)
            is_disease = any(kw in full_text for kw in disease_keywords)
            
            if is_control and not is_disease:
                control_samples.append(gsm_id)
            elif is_disease and not is_control:
                disease_samples.append(gsm_id)
            elif is_disease and is_control:
                # Ambiguous - skip
                pass
            else:
                # Default: try to split by source
                disease_samples.append(gsm_id)
        
        # If auto-detection fails, split samples 50/50
        if len(disease_samples) < 2 or len(control_samples) < 2:
            all_samples = list(gse.gsms.keys())
            mid = len(all_samples) // 2
            control_samples = all_samples[:mid]
            disease_samples = all_samples[mid:]
        
        return {
            "disease": disease_samples,
            "control": control_samples
        }

    def retrieve_datasets_for_disease(self, disease_name: str, n_datasets: int = 3) -> list:
        """
        Full pipeline to retrieve datasets for a disease.
        
        Args:
            disease_name: Disease name
            n_datasets: How many datasets to download
            
        Returns:
            List of processed dataset dicts
        """
        # Step 1: Search
        geo_ids = self.search_geo_datasets(disease_name, max_results=20)
        
        if not geo_ids:
            print("❌ No datasets found. Try a different disease name.")
            return []
        
        # Step 2: Get metadata
        datasets_meta = self.fetch_dataset_metadata(geo_ids)
        
        if not datasets_meta:
            print("❌ Could not retrieve dataset metadata.")
            return []
        
        print(f"\n📊 Top datasets found:")
        for i, d in enumerate(datasets_meta[:5], 1):
            print(f"  {i}. {d['gse_id']} - {d['n_samples']} samples - {d['title'][:60]}...")
        
        # Step 3: Download top datasets
        processed = []
        print(f"\n⬇️  Downloading top {n_datasets} datasets...")
        
        for meta in datasets_meta[:n_datasets + 3]:  # Try a few extra in case of failures
            if len(processed) >= n_datasets:
                break
                
            result = self.download_and_parse_gse(meta["gse_id"])
            
            if result and result.get("sample_groups"):
                groups = result["sample_groups"]
                if len(groups.get("disease", [])) >= 2 and len(groups.get("control", [])) >= 2:
                    processed.append(result)
                    print(f"  ✓ {meta['gse_id']}: {len(groups['disease'])} disease, {len(groups['control'])} control samples")
        
        print(f"\n✅ Successfully loaded {len(processed)} datasets")
        return processed
