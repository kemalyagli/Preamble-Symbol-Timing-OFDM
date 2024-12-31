# Preamble-Based Symbol Timing Algorithms in OFDM Systems

This repository contains the implementation and resources for the master's thesis project titled **"Preamble-Based Symbol Timing Algorithms in OFDM Systems"**. The project focuses on analyzing and comparing different preamble-based synchronization techniques in OFDM systems under various channel conditions.

---

## Repository Structure

```
Preamble-Based-Symbol-Timing-Algorithms-In-OFDM-Systems/
├── data/
│   ├── awgn_rayleigh_64_16.mat
│   ├── awgn_rayleigh_128_16.mat
│   └── awgn_rayleigh_256_32.mat
├── functions/
│   ├── metrics/
│   │   ├── schmidl_sync_metric.m
│   │   ├── minn_sync_metric.m
│   │   ├── park_sync_metric.m
│   │   ├── kim_sync_metric.m
│   │   └── ren_sync_metric.m
│   ├── preambles/
│       ├── schmidl_preamble.m
│       ├── minn_preamble.m
│       ├── park_preamble.m
│       ├── kim_preamble.m
│       ├── ren_preamble.m
│       └── proposed_preamble.m
├── main_scripts/
│   ├── main_awgn_channel.m
│   └── main_rayleigh_channel.m
├── publication/
│   └── publication.pdf
├── results/
│   └── results.m
└── thesis/
    ├── thesis.pdf
```

### 1. `data/`
Contains `.mat` files with simulated channel data for AWGN and Rayleigh fading channels, used in testing and analyzing synchronization algorithms:
- `awgn_rayleigh_64_16.mat`
- `awgn_rayleigh_128_16.mat`
- `awgn_rayleigh_256_32.mat`

### 2. `functions/`
#### `metrics/`
Functions implementing various synchronization metrics:
- **`schmidl_sync_metric.m`**: Implements the Schmidl's metric.
- **`minn_sync_metric.m`**: Implements the Minn's metric.
- **`park_sync_metric.m`**: Implements the Park's metric.
- **`kim_sync_metric.m`**: Implements the Kim's metric.
- **`ren_sync_metric.m`**: Implements the Ren's metric.

#### `preambles/`
Functions generating preambles for the synchronization algorithms:
- **`schmidl_preamble.m`**: Generates the Schmidl's preamble.
- **`minn_preamble.m`**: Generates the Minn's preamble.
- **`park_preamble.m`**: Generates the Park's preamble.
- **`kim_preamble.m`**: Generates the Kim's preamble.
- **`ren_preamble.m`**: Generates the Ren's preamble.
- **`proposed_preamble.m`**: Generates the proposed preamble introduced in the thesis.

### 3. `main_scripts/`
Scripts to test and analyze synchronization algorithms under different channel conditions:
- **`main_awgn_channel.m`**: Simulates and analyzes algorithms in an AWGN channel.
- **`main_rayleigh_channel.m`**: Simulates and analyzes algorithms in a Rayleigh fading channel.

### 4. `publication/`
Contains the PDF of the published paper:
- **`publication.pdf`**: The published research paper derived from the thesis.

### 5. `results/`
Contains the MATLAB script for post-processing and visualizing results:
- **`results.m`**: Processes simulation outputs and generates plots.

### 6. `thesis/`
Contains the PDF of the full thesis document:
- **`thesis.pdf`**: The complete thesis, including methodology, results, and conclusions.

---

## How to Use

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/Preamble-Based-Symbol-Timing-Algorithms-In-OFDM-Systems.git
   ```

2. Navigate to the repository:
   ```bash
   cd Preamble-Based-Symbol-Timing-Algorithms-In-OFDM-Systems
   ```

3. Open MATLAB and set the current directory to the repository folder.

4. Run the main scripts in the `main_scripts/` folder to simulate the algorithms:
   - `main_awgn_channel.m`
   - `main_rayleigh_channel.m`

5. Use `results/results.m` to process and visualize the results.

---

## Requirements

- MATLAB (R2020b or later recommended)
- Signal Processing Toolbox

---

## License
Feel free to use and modify the code for academic and research purposes. Please cite the publication if you use this work in your research.

---

## Citation
If you use this repository, please cite the publication:

> Yağlı, K., Aldırmaz Çolak, S., (2022). Preamble-Based Symbol Timing Algorithms in OFDM Systems. The European 
Journal of Research and Development, 2(2), 445– 458. DOI: 10.56038/ejrnd.v2i2.91

---

## Acknowledgments
This repository is part of the master's thesis project conducted at Kocaeli University. Special thanks to my advisor, Sultan Aldırmaz Çolak, and colleagues who provided guidance and support throughout this work.
