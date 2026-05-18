# Supplementary Code for TITLE

This repository contains the supplementary code used for the bioinformatic and imaging analyses described in the manuscript:

**TITLE**  
DOI: **XXX**

The study investigates chromatin organization, accelerated aging, transcriptomic alterations, and nuclear morphology in schizophrenia.

## Repository Structure

### `TranscriptomicClocks/`

Code for constructing transcriptomic aging clocks using Elastic Net regression.

**Data source:**
- PsychENCODE Consortium (DOI: 10.15154/1g4m-dy13)

---

### `AgeScore/`

Code for calculating the custom transcriptomic Age Score, a metric for estimating biological aging based on age-associated gene expression changes.

**Data sources:**
- Mendizabal et al., 2019 (DOI: 10.1186/s13059-019-1747-7)
- Zhu et al., 2024 (DOI: 10.7554/eLife.92393)

---

### `NucleiMicroscopy/`

Scripts for microscopy image processing and quantitative analysis of nuclear morphology.

**Data source:**
- Zenodo dataset (DOI: 10.5281/zenodo.20098460)

---

### `FineFeatureAnalysis/`

Code for analyzing fine-scale chromatin features: topologically associating domains (TADs) and chromatin loops in neuronal and non-neuronal cells from schizophrenia patients and healthy controls.

---

### `SNPContactAnalysis/`

Code for detection of significant, variant-sensitive, interactions and their downstream analysis.

---

## Requirements

The code is primarily written in:
- Python
- R
- Bash

Dependencies are specified within individual scripts and notebooks.

---

## Citation

If you use this code or the associated datasets, please cite:

**TITLE**  
DOI: **XXX**

---

## Contact

For questions regarding the code or datasets, please contact the corresponding authors listed in the manuscript.

---
## References

Akbarian, S., Liu, C., Knowles, J. A., Vaccarino, F. M., Farnham, P. J., Crawford, G. E., Jaffe, A. E., Pinto, D., Dracheva, S., Geschwind, D. H., Mill, J., Nairn, A. C., Abyzov, A., Pochareddy, S., Prabhakar, S., Weissman, S., Sullivan, P. F., State, M. W., Weng, Z., � Sestan, N. (2015). The PsychENCODE project. Nature Neuroscience, 18(12), 1707�1712. https://doi.org/10.1038/nn.4156 

Pletenev, I. A., Bazarevich, M., Zagirova, D. R., Kononkova, A. D., Cherkasov, A. V., Efimova, O. I., Tiukacheva, E. A., Morozov, K. V., Ulianov, K. A., Komkov, D., Tvorogova, A. V., Golimbet, V. E., Kondratyev, N. V., Razin, S. V., Khaitovich, P., Ulianov, S. V., & Khrameeva, E. E. (2024). Extensive long-range polycomb interactions and weak compartmentalization are hallmarks of human neuronal 3D genome. Nucleic Acids Research, 52(11), 6234�6252. https://doi.org/10.1093/nar/gkae271 

Pletenev, I. A., Vaulin, N., Molodova, M. N., Kuznechenkova, E., Soldatenkova, A., Efimova, O. I., Tvorogova, A. V., Khaitovich, P., Razin, S. V., Ulianov, S. V., & Khrameeva, E. E. (2025). Ultra-long-range Polycomb-coupled interactions underlie subtype identity of human cortical neurons (p. 2025.11.05.686502). bioRxiv. https://doi.org/10.1101/2025.11.05.686502 

Weng, Z. (2023). PsychENCODE Consortium. NIMH Data Repositories. https://doi.org/10.15154/1G4M-DY13 

Zagirova, D. R., Kononkova, A. D., Morozov, K. V., Molodova, M. N., Vaulin, N. S., Dudkovskaia, A. V., Dozorova, P. I., Efimova, O. I., Tvorogova, A. V., Ulianov, K. A., Khaitovich, P. E., Razin, S. V., Lagarkova, M. A., Ulianov, S. V., & Khrameeva, E. E. (2025). Fetal signatures in the 3D genome of iPSC-derived neurons: Implications for disease modeling (p. 2025.08.03.667702). bioRxiv. https://doi.org/10.1101/2025.08.03.667702 

Zhu, B., Ainsworth, R. I., Wang, Z., Liu, Z., Sierra, S., Deng, C., Callado, L. F., Meana, J. J., Wang, W., Lu, C., & Gonzalez-Maeso, J. (2024). Antipsychotic-induced epigenomic reorganization in frontal cortex of individuals with schizophrenia. eLife, 12, RP92393. https://doi.org/10.7554/eLife.92393
