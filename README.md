# D<sup>3</sup>E: Discrete Distributional Differential Expression

##About
D<sup>3</sup>E is a tool for identifying differentially-expressed genes, based on single-cell RNA-seq data. The main assumption we use in our method is that expression of transcripts follows Poisson-beta distribution with three parameters &alpha;, &beta; and &gamma;. The first two correspond to the rate of gene activation and deactivation, and the latter corresponds to the rate of transcription when a gene is in an active state.

##Files
- **D3EUtil.py** : a collection of methods for running the D<sup>3</sup>E analysis
- **D3ECmd.py** : an actual D<sup>3</sup>E tool (requires **D3EUtil.py** for execution)
- **D3EWed.py** : a light version of D<sup>3</sup>E, for running as a web-service (requires **D3EUtil.py** for execution)