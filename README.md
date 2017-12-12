# Fast Calculation for Galaxy Two-Point Correlations

A Python implementation of the algorithm described in [Paper Reference].

## Getting Started

Instructions to run the project on your local machine.

### Prerequisites

This implementation depends on the following Python packages:
* [AstroPy](http://www.astropy.org)
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

If desired, use virtualenv to install the prerequisites:
```
virtualenv venv
source venv/bin/activate
pip install numpy astropy scipy
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

We can quickly find the two-point correlation function xi for
galaxies in the DR9 survey of the souther sky using the included
configuration `configs/cmassS_coarse.py`.  This configuration differs
from `configs/cmassS.py` by using fewer bins and not employing any
strategies to reduce the number of evaluated galaxy pairs.  The
processing is divided into stages, each with their own corresponding
routine:
* `preprocessing.py` histograms the data from the galaxy catalogs
* `combinatorial.py` makes distributions of galaxy pairs
* `integration.py` incorporates cosmological data to convert angles and redshift to distances

### Preprocessing
Run the preprocessing routine, which takes seconds.
```
./lasspia.py configs/cmassS_coarse.py routines/preprocessing.py
```
View the headers of the output file.
```
./lasspia.py configs/cmassS_coarse.py routines/preprocessing.py --show
```
Run the preprocessing.plot() method to plot the contents of the output file (PDF output)
```
./lasspia.py configs/cmassS_coarse.py routines/preprocessing.py --plot
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
View the headers of the output file.
```
./lasspia.py configs/cmassS_coarse.py routines/combinatorial.py --show
```

### Integration
Run the integration routine, which takes seconds.
```
./lasspia.py configs/cmassS_coarse.py routines/integration.py
```
If you get "MemoryError", you can break the integration into slices of
bins of theta by passing `--nJobs` and `--nCores` (or `--iJob`)
arguments, and combining the output as in the prior step example.

View the headers of the output file.
```
./lasspia.py configs/cmassS_coarse.py routines/integration.py --show
```
Run the integration.plot() method.
```
./lasspia.py configs/cmassS_coarse.py routines/integration.py --plot
```

## Parallel and Batch Processing

The `combinatorial` and `integration` routines respond to the option
`--nJobs` to break processing into parallelizable units.  Combine with
the `--nCores` option to specify the number of processes to run
locally in parallel, or with the `--iJob` option to specify one (or
several) units to process.  Some batch systems which accept "job
arrays" use an environment variable to specify the job index, in which
case it may be convenient to specify `--iJobEnv` rather than `--iJob`.
After all units have succeeded, their outputs may be combined by
repeating the `--nJobs` command, absent any of `--nCores`, `--iJob`,
or `--iJobEnv`.  A minimal Dockerfile is included for use with batch
systems that run containers.

## Contributing

Pending

## Versioning

Pending

## Authors

* **Burton Betchart** - *Initial work* - [betchart](https://github.com/betchart)

See also the list of [contributors](https://github.com/betchart/lasspia/contributors) who participated in this project.

## License

Pending
