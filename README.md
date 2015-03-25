#D<sup>3</sup>E: Discrete Distributional Differential Expression

##About
D<sup>3</sup>E is a tool for identifying differentially-expressed genes, based on single-cell RNA-seq data. The main assumption we use in our method is that expression of transcripts follows Poisson-beta distribution with three parameters &alpha;, &beta; and &gamma;. The first two correspond to the rate of gene activation and deactivation, and the latter corresponds to the rate of transcription when a gene is in an active state. A web-version can be accessed at http://www.sanger.ac.uk/sanger/GeneRegulation_D3E/

##Files
- **D3EUtil.py** : a collection of methods for running the D<sup>3</sup>E analysis
- **D3ECmd.py** : an actual D<sup>3</sup>E tool (requires **D3EUtil.py** for execution)
- **D3EWed.py** : a light version of D<sup>3</sup>E, for running as a web-service (requires **D3EUtil.py** for execution)

##Running a Command Line Tool

To run D<sup>3</sup>E analysis please run D3ECmd.py:

```
python D3ECmd.py InputFile OutputFile Label1 Label2 [-m {1,2}] [-z {0,1}] [-n {1,0}] [-s {1,0}] [-v]
```

**Mandatory arguments:**

- InputFile : a path to the read-count table
- OutputFile : a path to the output file
- Label1 : a common label for the cells of the first type
- Label2 : a common label for the cells of the second type

**Optional arguments**

- -m : a run mode
- -z (--removeZeos) : if -z is set to 1, all zero enties in the read-count table will be ignored (default)
- -n (--normalise) : if -n is set to 1, a normalisation routine will be performed before analysis (default)
- -v : if -v flag is set, D3ECmd.py will run in a verbose mode

##Input File Format

D3ECmd.py accepts a tab-separated read-count table, where rows correspond to genes, and columns correspond to individual cells. The file should have a header row which has the following tab-separated format:

```
"GeneID	Label<sub>1</sub>	Label<sub>2</sub>	Label<sub>3</sub>	... "
```

where L<sub>i</sub> are the cell type labels. Differential expression analysis can be performed on two cell types at a time. 

Each line should start with a gene ID, followed by a sequence of read-counts. Empty lines are ignored.

##Run Modes

D3ECmd.py can run in three modes:

- Mode 1 : Methods of moments is used to estimate parameters ( fast but less accurate )
- Mode 2 : Bayesian inference is used to estimate parameters ( slow but more accurate, default )