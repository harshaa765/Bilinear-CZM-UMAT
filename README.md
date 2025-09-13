# 🔗 Bilinear Cohesive Zone Model (CZM) – Abaqus UMAT

This repository provides a **Fortran UMAT subroutine** for the **mixed-mode bilinear cohesive zone model (CZM)**, along with benchmark **Abaqus input files** for fracture simulations.  
It is intended as a reference and starting point for researchers and engineers working on cohesive fracture modeling in Abaqus.

---

## 📂 Repository Contents
- **`Bilinear_CZM_UMAT.for`** – Fortran subroutine implementing the mixed-mode bilinear cohesive zone model.  
- **Example Abaqus input files**:
  - `2D_DCB_Bilinear_CZM_UMAT.inp` – Double Cantilever Beam (Mode I).  
  - `2D_ENF_Bilinear_CZM_UMAT.inp` – End Notched Flexure (Mode II).  
  - `2D_MMB_Bilinear_CZM_UMAT.inp` – Mixed-Mode Bending specimen.  
  - `Job_1_Harsh_UMAT.inp` – Single element patch test.  

---

## ▶️ Usage
1. Compile the UMAT subroutine with Abaqus.  
2. Run any of the provided `.inp` files to reproduce standard fracture benchmarks.  

Example (4 CPUs, change paths as needed):
```bash
abaqus job=2D_DCB_Bilinear_CZM_UMAT user=Bilinear_CZM_UMAT.for cpus=4
```
---

## 🙌 Acknowledgment
If you find this repository useful, please ⭐ **star this repo** and **follow me** here on GitHub for more such content.  
Additionally, please cite the related research articles:

```bibtex
@article{sharma2025combined,
  title     = {Combined phase-field and cohesive zone modeling for mixed-mode fracture in polymer composites},
  author    = {Sharma, Harshdeep and Singh, Akhilendra},
  journal   = {Engineering with Computers},
  pages     = {1--29},
  year      = {2025},
  publisher = {Springer}
}

@article{sharma2024numerical,
  title     = {Numerical implementation of a modified cohesive zone model for HCF behavior of adhesively bonded composite laminates under mixed mode loading},
  author    = {Sharma, Harshdeep and Singh, Akhilendra},
  journal   = {International Journal of Fatigue},
  volume    = {181},
  pages     = {108128},
  year      = {2024},
  publisher = {Elsevier}
}
