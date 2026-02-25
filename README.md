# Canopy-SIP Model (Optical Version)

This repository contains the MATLAB implementation of the **Canopy-SIP Model** for simulating canopy optical Bidirectional Reflectance Factor (BRF).

## Overview

This model is designed for **discrete vegetation canopies**. It integrates:
1.  **Geometric-Optical (GO) Theory**: To calculate the area fractions of four scene components (sunlit/shaded crown, sunlit/shaded soil) and directional gap probabilities.
2.  **Spectral Invariants Theory (p-theory)**: To efficiently simulate multiple scattering within the canopy.


## Citation

If you use this code or model in your research, please cite the following paper which describes the model improvements and applications:
> **He, Y.**, Zeng, Y., Hao, D., Shabanov, N. V., Huang, J., Yin, G., Biriukova, K., Lu, W., Gao, Y., Celesti, M., Xu, B., Gao, S., Migliavacca, M., Li, J., & Rossini, M. (2025). Combining geometric-optical and spectral invariants theories for modeling canopy fluorescence anisotropy. *Remote Sensing of Environment*, 323, 114716. [https://doi.org/10.1016/j.rse.2025.114716](https://doi.org/10.1016/j.rse.2025.114716)

## Model Heritage & Acknowledgments

This model is built upon the foundational work of the original SIP model. We gratefully acknowledge the contributions of the original developers:

> **Zeng, Y.**, Xu, B., Yin, G., Wu, S., Hu, G., Yan, K., Yang, B., Song, W., & Li, J. (2018). Spectral Invariant Provides a Practical Modeling Approach for Future Biophysical Variable Estimations. *Remote Sensing*, 10(10), 1508. [https://doi.org/10.3390/rs10101508](https://doi.org/10.3390/rs10101508)

## Repository Structure

* `main.m`: The primary execution script. It sets up the sun-sensor geometry, loads structural parameters, runs the BRF simulation, and plots the results in the principal plane.
* `/CI_2/`: Directory containing structural look-up tables (gap fractions and clumping indices).
* `get_HSF_go.m`: Calculates the hotspot effect and bidirectional gap probabilities based on GO theory.
* `sunshade_H.m` / `sunshade_Kt_He.m`: Functions for calculating sunlit/shaded foliage and soil probabilities.
* `PHASE.m`: Calculates the scattering phase function.

## Usage

1.  Clone or download this repository.
2.  Open MATLAB and navigate to the repository folder.
3.  Run the `main.m` script.
4.  The script will:
    * Load necessary structural data.
    * Simulate BRF for a specific band (e.g., NIR) across different viewing angles.
    * Save the results to `BRF_SIP_SZA00_Nir2.mat`.
    * Automatically generate a plot showing the BRF distribution in the principal plane.

## License

This project is licensed under the MIT License.