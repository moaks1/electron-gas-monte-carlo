# Electron Gas Monte Carlo

A simple Monte Carlo electron-transport code for gas targets using LXCat-style collision cross-section data.

This project tracks electrons through a target gas or gas mixture. It samples collision times, chooses interaction channels from tabulated cross sections, updates particle energy and direction, and can create secondary electrons from ionization events.

The code is intentionally written in a basic, readable style so the physics and Monte Carlo logic are easy to follow.

---

## Features

- Electron transport through a gas target
- LXCat-style cross-section loading
- Elastic scattering
- Excitation energy loss
- Ionization with secondary-electron creation
- Attachment support, if attachment data are provided
- Pure gases and gas mixtures
- 3D particle histories
- Trajectory plotting from a Jupyter notebook
- Reproducible runs with a random seed

---

## Repository structure

A suggested layout is:

```text
electron-gas-monte-carlo/
├── README.md
├── requirements.txt
├── electron_mc_transport.py
├── electron_transport.ipynb
├── figures/
│   └── example_trajectory.png
└── data/
    └── README.md
```

The `data/` folder is not required to contain cross-section files in the public repository. See the note about data below.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR-USERNAME/electron-gas-monte-carlo.git
cd electron-gas-monte-carlo
```

Create and activate a conda environment:

```bash
conda create -n electron-mc python=3.11
conda activate electron-mc
```

Install the required packages:

```bash
pip install -r requirements.txt
```

Start JupyterLab:

```bash
jupyter lab
```

Then open:

```text
electron_transport.ipynb
```

---

## Cross-section data

This code expects LXCat-style electron collision cross-section data. A typical file contains blocks such as:

```text
ELASTIC
EXCITATION
IONIZATION
ATTACHMENT
```

with energy in eV and cross section in square meters.

Cross-section files are not included in this repository by default. Many LXCat datasets have their own citation and redistribution requirements, so it is better to download the data directly from the source and cite it properly.

For example, place your downloaded file in the repository folder and set:

```python
CROSS_SECTION_FILE = "eHexsec.txt"
```

inside the notebook.

---

## Basic usage

Load a single gas target:

```python
import numpy as np
import electron_mc_transport as emc

np.random.seed(1)

material = emc.load_lxcat_material(
    "eHexsec.txt",
    material_name="helium gas",
    target_name="He",
    number_density_m3=2.5e26,
)

electron = emc.ElectronParticle(
    energy_eV=100.0,
    material=material,
    idx=0,
)

secondary, event = electron.update(next_idx=1)

print(electron)
print(event)
```

---

## Gas mixture example

If cross-section data are available for each gas species, the code can also model a gas mixture.

```python
import electron_mc_transport as emc

n_total = 2.5e26

n2_processes = emc.load_lxcat_table("N2_xsecs.txt", target_name="N2")
o2_processes = emc.load_lxcat_table("O2_xsecs.txt", target_name="O2")
ar_processes = emc.load_lxcat_table("Ar_xsecs.txt", target_name="Ar")

air = emc.make_gas_mixture("air", [
    ["N2", 0.78 * n_total, n2_processes],
    ["O2", 0.21 * n_total, o2_processes],
    ["Ar", 0.01 * n_total, ar_processes],
])

electron = emc.ElectronParticle(
    energy_eV=100.0,
    material=air,
    idx=0,
)
```

For a mixture, the total macroscopic cross section is

```text
Sigma_total(E) = sum_i n_i sigma_i(E)
```

where `n_i` is the number density of species `i`, and `sigma_i(E)` is the microscopic cross section.

---

## Physics model

This is a direct, event-based Monte Carlo model. During each update:

1. The code evaluates the total macroscopic cross section at the current electron energy.
2. A collision time is sampled using a null-collision method.
3. The electron moves in a straight line during the sampled time step.
4. A collision process is sampled from the available cross sections.
5. The electron energy and direction are updated.
6. If ionization occurs, a secondary electron may be created.

The code currently assumes:

- no external electric or magnetic fields
- straight-line motion between collisions
- isotropic post-collision scattering
- nonrelativistic electron speed
- tabulated cross sections as the source of collision probabilities

---

## Limitations

This is a research/educational prototype, not a full production electron transport code.

Current limitations include:

- no electric-field acceleration
- no magnetic-field motion
- no differential angular cross sections
- no detailed secondary-electron energy spectrum
- no elastic energy loss to the target atom or molecule
- no condensed-matter or solid-state transport model
- no boundary geometry other than free particle histories

The code is best described as:

> Monte Carlo electron transport in gas targets using LXCat-style collision cross-section data.

It should not be described as a general-purpose solid-material electron transport code.

---

## Validation ideas

Good validation checks include:

1. Plot the loaded cross sections and compare them with the source data.
2. At fixed energy, verify that sampled collision types match the expected cross-section ratios.
3. Check that sampled free paths follow the expected exponential distribution.
4. Check that excitation removes the correct threshold energy.
5. Check that ionization conserves leftover kinetic energy after subtracting the ionization threshold.
6. For future electric-field versions, compare swarm quantities such as mean energy or drift velocity against Boltzmann or swarm-code results.

---

## Example output

The notebook can generate a 2D projection of the electron trajectory, showing the primary electron and any secondary electrons created through ionization.

Suggested figure location:

```text
figures/example_trajectory.png
```

---

## License

Choose a license before making the repository public. The MIT License is a common choice for a small educational code project.

Do not assume that downloaded cross-section data can be redistributed under the same license as your code. Keep the code license and the data citation/usage terms separate.

---

## Suggested citation

If you use LXCat or another database for electron collision cross sections, cite the database according to the source's recommended citation format.

You can also add a `CITATION.cff` file later if you want people to cite this GitHub repository directly.
