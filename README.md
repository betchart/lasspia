# Fast Calculation for Galaxy Two-Point Correlations

A Python implementation of the algorithm described in [Paper Reference].

## Getting Started

Instructions to run the project on your local machine.

### Prerequisites

This implementation depends on the following Python packages:
* [AstroPy (v1)](http://www.astropy.org)
* [SciPy library](https://github.com/scipy/scipy)
* [NumPy](http://www.numpy.org)
* [Matplotlib](http://matplotlib.org) (optional)

The included configurations rely on DR9 data from the Sloan Digital Sky Survey:
* [SDSS archive](https://data.sdss.org/sas/dr9/boss/lss/)

### Installing

Clone the repository to your local system:

```
git clone https://github.com/betchart/lasspia.git
```

Download the SDSS galaxy survey catalogs:

```
cd lasspia/
./downloadDR9.py
```

Test your installation by running the 'quickscan.py' routine.

```
./lasspia.py configs/cmassS.py routines/quickscan.py
```

## A Complete Example

We can quickly find the two-point correlation function xsi for
galaxies in the DR9 survey of the souther sky using the included
configuration `configs/cmassS_coarse.py`.  This configuration differs
from `configs/cmassS.py` by using fewer bins and not employing any
strategies to reduce the number of evaluated galaxy pairs.  The
processing is divided into stages, each with their own corresponding
routine; * 'preprocessing' histograms the data from the galaxy
catalogs * 'combinatorial' makes distributions of galaxy pairs *
'integration' incorporates cosmological data to convert angles and
redshift to distances

### Preprocessing
Run the preprocessing routine, which takes seconds.
```
./lasspia.py configs/cmassS_coarse.py routines/preprocessing.py
```
View the headers of the preprocessing output file.
```
./lasspia.py configs/cmassS_coarse.py routines/showPre.py
```

### Combinatorial
Run the combinatorial routine, which takes about half an hour.
```
./lasspia.py configs/cmassS_coarse.py routines/combinatorial.py
```
Alternatively, run the combinatorial routine with parallel jobs.
```
./lasspia.py configs/cmassS_coarse.py routines/combinatorial.py --nJobs 8 --nCores 4
```
Combine the output of the jobs.
```
./lasspia.py configs/cmassS_coarse.py routines/combinatorial.py --nJobs 8
```

### Integration
Run the integration routine, which takes seconds.
```
./lasspia.py configs/cmassS_coarse.py routines/integration.py
```

### Visualization
Graph the results
```
pending
```

## Contributing

Pending

## Versioning

Pending

## Authors

* **Burton Betchart** - *Initial work* - [betchart](https://github.com/betchart)

See also the list of [contributors](https://github.com/betchart/lasspia/contributors) who participated in this project.

## License

Pending
