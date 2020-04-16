# actsims

This package lets you generate simulations of ACT and Planck single frequency maps. The simulations are used for verification and covariance calculations in the ACT power spectrum, lensing and component separation pipelines.

Details of the simulations will be available in Choi. et. al. (in prep). The broad model is:
* lensed CMB on the full-sky
* Gaussian foregrounds on the full-sky
* beam convolution on the full-sky
* Gaussian noise model that has 2D Fourier space anisotropy (as an average across the patch) and spatial inhomogenity, but no inhomogenity of anisotropy

## Dependencies

As a user, you can generate simulations if you have:

* Code: pixell, soapack
* Data: lensed CMB harmonics data, noise template data


As a developer, you can generate new lensed CMB realizations and create noise templates if you have:

* Code: pixell, soapack, tilec

## Installation

To install: `pip install -e . --user`

## Usage (``user'' mode)


## Usage (``developer'' mode)


## Developers

Please contact Alex van Engelen, Dongwon 'DW' Han and Mat Madhavacheril for questions regarding the simulations.