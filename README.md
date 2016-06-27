# A Bioconductor workflow for processing and analysing spatial proteomics data

[![Build Status](https://travis-ci.org/lgatto/bioc-pRoloc-hyperLOPIT-workflow.svg?branch=master)](https://travis-ci.org/lgatto/bioc-pRoloc-hyperLOPIT-workflow)

NB: think about creating new `MSnSet` for the `hyperLOPIT2015` or writing a function to strip it down some of the `fData` for the purpose of demonstrating this workflow

### Data import
- Start with protein level and `readMSnSet2`
- Brief explanation of what a `MSnSet` is, how to access and what can be found in the `pData`, `fData`, refer to `hyperLOPIT2015` specifically (e.g. two replicates combined refer to MWBT paper, acquired using MS3, intrument etc.)

### Normalisation
- Normalise by sum
- Highlight that many methods available to use directly with `MSnSet` in `MSnbase`

### Quality control
- Data structure, quality control
- Mention `pRolocGUI` here

### Markers
- importance of organelle markers and implications
- `addMarkers`, `getMarkers` and friends
- manual curation, quality control
- revisting visualisation with markers

### Novelty detection

### Supervised ML
- thresholding

### Inductive transfer 

### Post-processing
- overlap with GO, enrichment study?
