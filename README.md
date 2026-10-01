# AMP Atlas — Vinay's top 100

**Live explorer:** https://aarnavbhat23.github.io/amp-top100-atlas/

Independent visualization by Aarnav Bhat, based on Vinay Prabhu's submission. This repository is separate from the original submission and does not modify it. GitHub Pages serves the static files from the root of the main branch.

This portable explorer contains all 100 Boltz-2.1 predicted monomer structures, sequence descriptors, reference comparisons, residue confidence and PAE matrices. It is an ESM Atlas–inspired local explorer, not affiliated with Meta's ESM Atlas.

## Open

Unzip the archive, open a terminal in this folder and run:

```sh
python3 -m http.server 8903 --bind 127.0.0.1
```

Open http://127.0.0.1:8903 in a WebGL-enabled browser. Opening index.html directly is insufficient because browsers restrict local-file data loading. The app's assets are bundled; it does not send sequences to a remote visualization service. Source links open external sites only when clicked.

## Explore

- Drag the map to rotate; scroll to zoom. Click a ribbon, choose a rank, or search a sequence substring.
- Change map colour to compare confidence, length, charge, rank or reference similarity.
- The right-hand molecular viewer rotates independently. Expand it, switch representations, or select individual sequence letters to highlight atoms.
- Scroll the detail panel for calculated properties, reference matches and PAE. Download original CIF, derived PDB or raw confidence metrics for each peptide.
- Methods & evidence explains the models, projection diagnostics and limitations.

## Evidence and limitations

Source commit: f77bfdbbf2aaaa4442d7453591d8a1c6ee224588 of vinayprabhu/amp-challenge-2027-submission. Exact CSV input SHA-256: 3cc4a8a7652213f9bd75073737601b924bb1745f68274e20dca966c4336a1b5b.

All 100 structures pass exact sequence, chain, coordinate and confidence validation. See structure_validation.json. PDB B factors represent predicted confidence, not experimental displacement factors. Native CIF files and raw metrics are retained in structures/.

Map: ESM-2 35M residue-mean embeddings → 3D cosine UMAP, seed 42. All settings and projection diagnostics are in atlas_method.json. Map spacing has no physical units and does not validate functional clusters. Only the top100 are embedded, not the complete library.

The reference ratio is the audited validator's normalized Levenshtein ratio, not alignment percent identity. Top40 has ratio 0.8000, allowed by the pinned strict >0.8 cutoff. Top60 has one cysteine; no disulfide was requested. Consult audit.json and validator_checks.json for the full sequence audit. Passing these checks is not complete submission or chemical eligibility certification.

These are predictions, not experimental structures or antimicrobial assays. Short peptides may change conformation in different environments. Molecular properties assume unmodified free termini. The app does not reproduce or validate the submission's training/ranking pipeline.

Rendering uses 3Dmol.js 2.5.3 (MIT; https://3dmol.org). Embeddings use facebook/esm2_t12_35M_UR50D. Raw predictions were obtained using approved paid Boltz jobs; the quoted combined estimate was $2.50, not verified billing.
