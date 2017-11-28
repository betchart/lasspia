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
git clone https://github.com/betchart/baofast.git
```

Download the SDSS galaxy survey catalogs:

```
cd baofast/
./downloadDR9.py
```

Test your installation by running the 'quickscan.py' routine.

```
./baofast.py configs/cmassS.py routines/quickscan.py
```

## A Complete Example

Pending

## Contributing

Pending

## Versioning

Pending

## Authors

* **Burton Betchart** - *Initial work* - [betchart](https://github.com/betchart)

See also the list of [contributors](https://github.com/betchart/baofast/contributors) who participated in this project.

## License

Pending
