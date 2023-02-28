<p align="center">
	<img src="img/logo.png" width="564" height="300">
</p>

# CUSTOM: Codon Usage to Specific Tissue OptiMizer

This package provides a codon optimization tool for tissue-specific gene design. It follows a probabilistic approach with two steps:
1) Translate tissue-specific codon preferences into a pool of optimal sequences
2) Select the desired sequence based on parameters of relevance.

## Dependencies

- Python >= 3.7
- [pandas](https://pandas.pydata.org/), [numpy](https://numpy.org/), [RNA](https://github.com/ViennaRNA/ViennaRNA)

## Installing

1) (optional) Create a virtual environment to install the tool

2) Install CUSTOM and its requirements using pip:
```bash
pip install custom_optimizer
```

## Basic usage

As a basic example, here is the code to optimize an eGFP protein to kidney:
```python
# Import package
import custom
# Start the optimizer
opt = TissueOptimizer("kidney", n_pool=100)
# Optimize the eGFP sequence
egfp = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
opt.optimize(egfp)
# Select the top 10 sequences
best_egfp_lung = opt.select_best(by={"MFE":"min","MFEini":"max","CAI":"max","CPB":"max","ENC":"min"},homopolymers=7, top=10)
```

## Contact

This project has been developed at [Center for Genomic Regulation](http://www.crg.eu/).

[Xavier Hernandez-Alias](mailto:xavier.hernandez@crg.eu)

[Martin H. Schaefer](mailto:martin.schaefer@ieo.it)

[Luis Serrano](mailto:luis.serrano@crg.eu)

[Schaefer laboratory](https://www.schaeferlab.org)

[Serrano laboratory](http://serranolab.crg.eu)


## Cite

Hernandez-Alias, X., Benisty, H., Radusky, L.G., Serrano, L. & Schaefer, M. H. (2023). Using protein-per-mRNA differences among human tissues in codon optimization. Genome Biology, 24(1):34. (https://doi.org/10.1186/s13059-023-02868-2)


## License

CUSTOM is under a common GNU GENERAL PUBLIC LICENSE. Plese, check [LICENSE](./LICENSE) for further information.

###### [2021] - Centre de Regulació Genòmica (CRG) - All Rights Reserved*
