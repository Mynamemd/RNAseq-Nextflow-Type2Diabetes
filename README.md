# 🧬 End-to-End Automated RNA-Seq Pipeline & Differential Gene Expression Analysis
*(Type 2 Diabetes Mellitus vs. Control Case Study)*

![Nextflow](https://img.shields.io/badge/Nextflow-26.04.6-brightgreen?style=for-the-badge&logo=nextflow)
![R](https://img.shields.io/badge/R-4.4+-blue?style=for-the-badge&logo=r)
![Bioconductor](https://img.shields.io/badge/Bioconductor-DESeq2-green?style=for-the-badge)
![Status](https://img.shields.io/badge/Pipeline-Production--Ready-orange?style=for-the-badge)

---

## 📌 Executive Summary

This repository contains a robust, scalable, and fully automated **RNA-Seq data processing pipeline** built using **Nextflow (DSL2)**, coupled with downstream **Differential Expression (DGE)** and **Functional Pathway Enrichment Analysis** in **R/Bioconductor**.

The pipeline was executed on paired-end transcriptomic datasets comparing **Type 2 Diabetes (T2D) patients against healthy controls**, successfully identifying key biological pathways, dysregulated genes, and mitochondrial metabolic disruptions associated with T2D progression.

---

## 🛠️ Pipeline Architecture & Workflow

The workflow integrates industry-standard bioinformatics utilities in a modular Nextflow DAG (Directed Acyclic Graph) architecture with checkpoint caching (`-resume` enabled):

```text
       ┌────────────────────────┐
       │   Raw Paired-End FASTQ │
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │     FastQC (QC)        │
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │     fastp (Trimming)   │
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │ HISAT2 (Genome Align)  │ ◄─── GRCh38 Genome Index
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │ featureCounts (Matrix) │ ◄─── Ensembl GTF Annotation
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │ DESeq2 & Enrichment    │ ───► Volcano Plots, GO & KEGG Pathways
       └────────────────────────┘
