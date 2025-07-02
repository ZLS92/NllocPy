# NllocPy

**NllocPy** is a Python wrapper for **NonLinLoc**, an open-source software developed by Anthony Lomax for probabilistic earthquake location.  
This package provides Python tools to automate input file preparation, manage control parameters, run NonLinLoc, and analyze or visualize results.

---

## 🚧 Project Status

**This project is currently under development.**  
Features, interfaces, and file structures may change. Use it at your own risk and feel free to contribute!


---

## ⚙️ Features

- **Input File Automation**: Easily create control parameter files (`nlloc_control_param.txt`), velocity models, and travel-time grids.
- **Utility Scripts**: Python tools for handling raster DEMs and geospatial shapefiles.
- **Integrated NonLinLoc Core**: Includes the **original NonLinLoc C source code** (`src/`) for easy compilation.
- **Result Visualization**: Python functions for plotting and interactive analysis.

---

## 👤 Attribution

- **NllocPy** was developed and is maintained by **Dr. Luigi Sante Zampa**.
- The **C code in `src/`** is the original **NonLinLoc** software developed by Anthony Lomax.  
  NonLinLoc is available at [ut-beg-texnet/NonLinLoc](https://github.com/ut-beg-texnet/NonLinLoc).
- **NllocPy** only provides a Python automation layer and utility scripts for easier integration.
- Some ideas and components for reading/writing NonLinLoc grid files uses the [NLLGrid](https://github.com/claudiodsf/nllgrid),  
  a Python library by Claudio Satriano and Natalia Poiata.
- If you use this wrapper, please **also cite NonLinLoc** in your scientific work:

  > Lomax, A., Virieux, J., Volant, P., & Berge, C. (2000).  
  > *Probabilistic earthquake location in 3D and layered models: Introduction of a Metropolis-Gibbs method and comparison with linear locations.*  
  > In **Advances in Seismic Event Location**, Springer, pp. 101–134.  
  > DOI: [10.1007/978-94-015-9536-0_5](https://doi.org/10.1007/978-94-015-9536-0_5)