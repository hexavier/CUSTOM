<p align="center">
	<img src="img/logo.png" width="371" height="300">
</p>

# TisOpt: Codon optimizer for tissue-specific expression

This package provides a codon optimization tool for tissue-specific gene design. It follows a probabilistic approach with two steps:
1) Translate tissue-specific codon preferences into a pool of optimal sequences
2) Select the desired sequence based on parameters of relevance.

## Dependencies

- Python >= 3.7
- [pandas](https://pandas.pydata.org/), [numpy](https://numpy.org/), [RNA](https://github.com/ViennaRNA/ViennaRNA)

## Installing

1) (optional) Create a virtual environment to install the tool

2) Install TisOpt and its requirements using pip:
```bash
pip install tisopt
```

## Basic usage

As a basic example, here is the code to optimize an eGFP protein to kidney:
```python
# Import package
import tisopt
# Start de optimizer
opt = TissueOptimizer("kidney", n_pool=100)
# Optimize the eGFP sequence
egfp = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
opt.optimize(egfp)
# Select the top 10 sequences
best_egfp_lung = opt.select_best(by={"MFE":"min","MFEini":"max","CAI":"max","CPB":"max","ENC":"min"},homopolymers=7, top=10)
```

## Contact

Xavier Hernandez-Alias: xavier.hernandez@crg.eu

Martin H. Schaefer: martin.schaefer@ieo.it

Schaefer laboratory: https://www.schaeferlab.org

Serrano laboratory: http://serranolab.crg.eu


## Cite

Citation text (https://doi.org/)
