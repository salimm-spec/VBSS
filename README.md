%% Khalid W. Hameed  && Salim A. Mohammed ALi
MATLAB codes run on Matlab version 2021a

# Physical Layer Security Simulation: VBSS vs. FHSS under Blind Eavesdropping

A MATLAB simulation framework evaluating the physical-layer security, bit error rate (BER), and theoretical secrecy rate of **Variable Bandwidth Subband Spreading (VBSS)** compared against traditional **Frequency Hopping Spread Spectrum (FHSS)** over OFDM/Rayleigh fading channels.

The framework models legitimate multi-user communication alongside a passive, blind eavesdropper (Eve) utilizing energy detection and unsupervised machine learning (K-means clustering) to intercept, associate, and demodulate user transmissions without prior knowledge of the hopping/allocation tables.

---

## Repository Structure

| File | Type | Description |
| :--- | :--- | :--- |
| `VBSS_Observation_V8a.m` | Main Script | Master Monte Carlo simulation comparing VBSS and FHSS across varying user counts and SNR levels, generating BER and eavesdropper attack plots. |
| `VBSS_Table3.m` | Function | Generates dynamic subcarrier allocation tables (subband widths `q` and starting positions `s`) across hops while guaranteeing target average bandwidth `qavg`. |
| `EavesdropObservation.m` | Function | Simulates a blind eavesdropper performing MAD-based energy detection, log-feature extraction, K-means clustering, optimal permutation user-matching, and 16-QAM demodulation. |
| `Chan_Ray_AWGN.m` | Function | Passes time-domain signals through a 7-tap frequency-selective Rayleigh fading channel with Additive White Gaussian Noise (AWGN). |
| `RX_module.m` | Function | Performs cyclic prefix (CP) removal, FFT conversion to the frequency domain, and Zero-Forcing (ZF) channel equalization. |
| `conf.m` | Function | Computes two-sided Student's $t$-distribution confidence intervals across Monte Carlo trials. |
| `RS_curves_2.m` | Standalone Script | Analytically computes and plots theoretical secrecy rates ($R_s$) across varying subcarrier counts ($n$), user counts ($N$), and SNRs. |

---

## Prerequisites & Toolboxes

* **MATLAB** (R2021a )
* **Communications Toolbox** (`comm.RectangularQAMModulator`, `comm.RectangularQAMDemodulator`, `awgn`)
* **Statistics and Machine Learning Toolbox** (`kmeans`, `tinv`, `statset`)

---

## How It Works

### 1. Resource Allocation (`VBSS_Table3.m`)
* **VBSS Mode:** Dynamically varies each user's subband width per hop within $[q_{\min}, q_{\max}]$ while enforcing an exact cumulative carrier budget ($C_{\text{target}} = q_{\text{avg}} \times T$) across $T$ hops, randomized with non-overlapping frequency gaps.
* **FHSS Mode:** Reuses `VBSS_Table3` with fixed bandwidth parameters ($q_{\min} = q_{\max} = q_{\text{avg}} = 10$), yielding constant-width frequency hopping.

### 2. OFDM Transceiver & Channel (`VBSS_Observation_V8a.m`, `Chan_Ray_AWGN.m`, `RX_module.m`)
* User bitstreams are modulated using **16-QAM** (`4 bits/subcarrier`), mapped onto allocated subbands inside a 512-point IFFT, and prepended with a 20-sample cyclic prefix.
* Signals traverse a normalized 7-tap Rayleigh multipath channel (`[1 0.7 0.5 0.35 0.25 0.15 0.1]`) plus AWGN.
* Legitimate receivers equalize the channel using perfect channel state information (CSI) and extract their assigned subbands.

### 3. Blind Eavesdropper Model (`EavesdropObservation.m`)
Eve receives the composite OFDM signal without access to the hopping/allocation tables and executes a four-stage attack:
1. **Blind Energy Detection:** Estimates the noise floor via Median Absolute Deviation (MAD) on smoothed power spectra ($P_{\text{sm}}$) to isolate active frequency bands.
2. **Feature Extraction & Clustering:** Extracts three observable features per detected band—$\log(1 + \text{width})$, $\log(1 + \text{energy})$, and $\log(1 + \text{energy}/\text{width})$—standardizes them, and groups observations into $N_{\text{users}}$ anonymous clusters via K-means (50 replicates).
3. **Evaluation Association:** Matches detected bands to ground-truth allocations using an Intersection-over-Union ($\text{IoU} \ge 0.50$) threshold and solves the optimal cluster-to-user permutation.
4. **Message Recovery:** Demodulates 16-QAM symbols from correctly associated bands to evaluate Eve's BER and overall message recovery rate.

---

## Usage

### Running the Full Link-Level Simulation
1. Open MATLAB and navigate to the repository directory.
2. Open `VBSS_Observation_V8a.m` and adjust simulation parameters at the top of the script if desired:
   * `UserCount = [2 3 5 8];` — Number of simultaneous users
   * `snr_range = -5:10;` — SNR range in dB
   * `xx = [50];` — Number of Monte Carlo runs (`MC`)
   * `N_hops = 200;` — Number of frequency hops per run
3. Run the script:
   ```matlab
   VBSS_Observation_V8a
