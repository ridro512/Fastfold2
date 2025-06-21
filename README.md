# FastFold 2 - Streamlined MSA Generation for AlphaFold Users

FastFold 2 is a lightweight modification of the AlphaFold 2.3.1 (non-docker) pipeline. It replaces the resource-intensive HHblits step with MMseqs2 easy-search, significantly improving speed and usability for high-performance computing environments.

This fork was made  for users running AlphaFold on shared HPC systems like Australia's Gadi (NCI).

Made by K. Zhao.

# 🚧 Disclaimer
This is an experimental modification. It has been tested on Gadi (NCI HPC), but your mileage may vary depending on system configuration and database setup.

# 🔐 License
This project is based on code licensed under the Apache License 2.0. See LICENSE for full terms.

AlphaFold model parameters are released under the CC BY 4.0 license.

# 📘 Background
AlphaFold is a breakthrough in protein structure prediction developed by DeepMind. This modified implementation does not change the model weights or inference code. Instead, it offers an alternative MSA generation pipeline for users who prefer MMseqs2 over HH-suite tools.

For academic use, please continue to cite the original works:

Jumper et al., Nature (2021), https://doi.org/10.1038/s41586-021-03819-2

Evans et al., bioRxiv (2021), https://doi.org/10.1101/2021.10.04.463034
