# University of Michigan Bioinformatics Core
![DNA](res_brcf_bioinformatics_dna_stock_blue.jpeg)

We help researchers identify and interpret patterns in RNA, DNA, and proteins by placing high-throuhgput data into a biologically meaningful context. This includes:

- **Elaborating experimental design**
- **Developing reproducible workflows**
- **Analyzing next-generation sequencing data**
- **Assisting in visualizing and interpreting results**
- **Letters of support, grant proposal development, or manuscript revision**
- **1-1 project analysis support and training**
- **Virtual, hands-on bioinformatics workshops** (see [current workshop schedule and descriptions](https://michmed.org/XYQwq))
- **Recharge and percent-effort support models**

## We would love to talk with you

- See [our website](https://michmed.org/GqGzZ) for more information about the core staff.
- Set up a [free consulation](https://docs.google.com/forms/d/e/1FAIpQLSepk7VqOl3xmBgkZybrl71VuQmKk3YmkgmpaBO4dD2hOtIh4w/viewform) to discuss your idea/project.
- Send us an email at bioinformatics@umich.edu .

## Bioinformatics analyses we support

<!-- 
Note:
- The <details> id attribute below provides a bookmark to expand the section from a URL. 
- Github will prefix this attribute with "user-content-" to avoid potential collisions 
  with Markdown generated HTML.
- The <summary> line is linked to the preceding details@id (with the added prefix) which 
  makes it easy link to the expanded section.
- Making the <summaries> links instead of headings also improves the spacing and makes the 
  expand/collapse indicator responsive (i.e. you can click on the indicator and it will 
  actually expand/collapse).
- Many Bothans died to bring us this information.
-->

<details id='expand-bulk-rna-seq-gene-expression'>
<summary><a href='#user-content-expand-bulk-rna-seq-gene-expression'>Bulk RNA-Seq gene expression</a></summary>

- Poly(A) selection, total RNA, miRNA, Ribo-Seq, long-read gene expression
- Differential gene expression
- Differential isoform expression, isoform switching
- Allele specific expression 
- Functional enrichment analysis (GO terms, KEGG pathways)
- [Sample RNA-Seq analysis report](https://umich-brcf-bioinf.github.io/Watermelon/doc/SampleReport.html)
- Tools & Resources:

   - [Workshop: RNA-Seq Demystified](https://medresearch.umich.edu/office-research/about-office-research/biomedical-research-core-facilities/bioinformatics-core/bioinformatics-workshops-training#rna-seq-demystified)
   - [nf-core/rnaseq analysis pipeline](https://nf-co.re/rnaseq) | [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) | [iPathwayGuide](https://advaitabio.com/bioinformatics/ipathwayguide/) | [WebGestalt](https://www.webgestalt.org/) | [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp)
<hr/>
</details>

<details id='expand-single-cellsingle-nuclei-analysis'>
<summary><a href='#user-content-expand-single-cellsingle-nuclei-analysis'>Single-cell/single-nuclei analysis</a></summary>
 
- scRNA-Seq/ snRNA-Seq gene expression: 3', 5', Flex
- V(D)J immune profiling
- Cell surface protein profiling, CITE-seq (a.k.a. TotalSeq, ADT)
- snATAC-Seq
- snRNA-Seq + snATAC-Seq
- Trajectory analysis, Velocity analysis
- Single-cell analysis of long-reads
- Tools & Resources:
  
  - [Workshop: Intro to Single-Cell Analysis](https://medresearch.umich.edu/office-research/about-office-research/biomedical-research-core-facilities/bioinformatics-core/bioinformatics-workshops-training#intro-to-single-cell-analysis)
  - [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest)] | [Seurat](https://satijalab.org/seurat/) | [scCatch](https://github.com/ZJUFanLab/scCATCH) | [Monocle](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/) | [veloctyo](https://velocyto.org/)
<hr/>
</details>

<details id='expand-spatial-sequencing-analysis'>
<summary><a href='#user-content-expand-spatial-sequencing-analysis'>Spatial sequencing analysis</a></summary>
 
- Visium / Visium HD
- Xenium in-situ/ subcellular
- GeoMX DSP
- Tools & Resources:

  - [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger/latest) | [Seurat](https://satijalab.org/seurat/) | [Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer/latest) | [GeoMX tools](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html)
<hr/>
</details>

<details id='expand-epigenomic-analysis'>
<summary><a href='#user-content-expand-epigenomic-analysis'>Epigenomic analysis</a></summary>
 
- DNA Methylation from WGBS/oxBS/EM-Seq, ERRBS/oxERRBS, long-reads
- Chromatin accessibility from bulk ATAC-Seq
- Histone profiling from ChIP-Seq / Cut & Run / Cut & Tag
- Transcription factor binding from ChIP-Seq
- EPIC-Array
- Chromatin conformation from Hi-C, 3C, ChIA-PET
- biomodal evoC

- Tools & Resources:

  - [nf-core/methylseq](https://nf-co.re/methylseq) | [nf-core/atacseq](https://nf-co.re/atacseq) | [nf-core/chipseq](https://nf-co.re/chipseq) | [nf-core/cutandrun](https://nf-co.re/cutandrun) | [HiCUP](https://www.bioinformatics.babraham.ac.uk/projects/hicup/) | [Juicer](https://github.com/theaidenlab/juicer) | [HiC-Pro](https://github.com/nservant/HiC-Pro
) | [Biomodal Duet](https://biomodal.com/technology/duet-software/) | [annotatr](https://bioconductor.org/packages/devel/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html) | [methylSig](https://www.bioconductor.org/packages/release/bioc/html/methylSig.html)
<hr/>
</details>

<details id='expand-variant-identification-assembly'>
<summary><a href='#user-content-expand-variant-identification-assembly'>Variant identification / assembly</a></summary>
 
- Variant identification / structural variation from WGS, exome, panel, long-reads
- Germline & somatic variants
- Copy Number Analysis from WGS
- Variant impact annotation
- Genome assembly from short reads, long-reads, hybrid. 
- Transcriptome assembly
- Tools & Resources:

  - [nf-core/sarek](https://nf-co.re/sarek) | [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) | [SnpEff](https://pcingola.github.io/SnpEff) | [SPAdes](https://github.com/ablab/spades) | [Velvet](https://github.com/dzerbino/velvet) | [Flye](https://github.com/mikolmogorov/Flye) 
<hr/>
</details>

<details id='expand-microbiome-analysis'>
<summary><a href='#user-content-expand-microbiome-analysis'>Microbiome analysis</a></summary>
 
- 16S amplicon
- Metagenomics / Metatranscriptomics
- Transposon sequencing (Tn-seq)
- Tools & Resources:

  - [mothur](https://mothur.org/) | [SqueezeMeta](https://github.com/jtamames/SqueezeMeta)

<hr/>
</details>

<details id='expand-proteome-expression-analysis'>
<summary><a href='#user-content-expand-proteome-expression-analysis'>Proteome expression analysis</a></summary>
 
- Liquid Chromatography-Mass Spectrometry (label free, isobaric labeling/TMT, DIA)
- Antibody/Aptamer (ELISA, OLink, Somalogic)
- Tools & Resources:

   - [Perseus](https://www.maxquant.org/perseus/) | [MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html) | [STRING](https://string-db.org/) | | [iPathwayGuide](https://advaitabio.com/bioinformatics/ipathwayguide/) | [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) | [Olink Insight](https://olink.com/software/olink-insight) | [Olink Analyze](https://github.com/Olink-Proteomics/OlinkRPackage)
<hr/>
</details>
