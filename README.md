# tempo_obs_SiFAP2.py

A lightweight Python script to compute the required observation exposure time for **SiFAP2** (Silicon Fast Optical Photometer 2) proposals at the Telescopio Nazionale Galileo (TNG).

The script calculates the total observation time ($T_{\text{obs}}$) needed to detect a specific background-subtracted pulse amplitude at a given significance level, incorporating calibration data, statistical frameworks from Vaughan et al. (1995), and corrections for current detector crosstalk levels.

---

## Features

* **Calibration Integration:** Automatically estimates net source count rates ($R_{\text{net}}$) based on the official SiFAP2 V-band magnitude calibration curve.
* **Statistical Robustness:** Uses the Vaughan et al. (1995) framework to handle single-trial or multi-trial period searches.
* **Crosstalk Correction:** Includes an updated power threshold correction to account for current detector crosstalk levels (adjusting the single-trial power threshold from $\sim11.83$ to $\sim13.25$ at a 0.0027 false-alarm probability).
* **Double-Check Verification:** Recalculates the expected background-subtracted amplitude using the derived $T_{\text{obs}}$ to ensure mathematical consistency.

---

## Prerequisites & Installation

### Requirements
* Python 3.x
* Standard Python libraries (no heavy external dependencies required).

### Setup
1. Download `tempo_obs_SiFAP2.py` and place it into your working directory.
2. Make the script executable via your terminal:
   ```bash
   chmod +x tempo_obs_SiFAP2.py
