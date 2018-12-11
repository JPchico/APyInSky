# **A**tomistic **Py**thon **In**tegrated **Sk**yrmion analysis framework (APyInSky)

The `APyInSky` package is a set of `python3` scripts aimed to the determination of several relevant quantities for skyrmion dynamics. It aims to bring under one umbrella the determination of several key quantities of magnetic skyrmions in a consistent fashion.

The scripts are designed to work with the output produced by the [`UppASD`](https://github.com/UppASD/UppASD "UppASD in Github") software package. Future support is planned for `ovf` files.

## Skyrmion dynamics analysis

`APyInSky`, currently is capable of calculating the following quantities:

- Skyrmion position as a function of time.
- Skyrmion velocity as a function of time.
    * From raw data.
    * From interpolated and smoothed data.
- Skyrmion radius as a function of time.
    * From raw data.
    * From interpolated and smoothed data. 
- Skyrmion circularity as a function of time.
- Skyrmion profile as a function of time.
- Determination of the magnetic forces acting over the skyrmion as a function of time.
- Topological charges.
- Determination of the dissipation matrix.

### Use
The control options for the package are specified in a file dubbed `Skx_inp.yml`. A minimum set of inputs, such as which kind of analysis must be done, path to the data, among others must be provided. 
To run the script is is enough to execute the main `APyInSky.py` script, found in the `src` folder.

In the `tools` folder a script to generate a basic input file is also present.

## Attributions
This package is developed by Jonathan Chico, with contributions from Imara Lima Fernandes and Markus Hoffmann.