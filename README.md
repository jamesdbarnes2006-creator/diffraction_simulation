# Diffraction Simulation

This project simulates diffraction through an aperture using scipy's dblquad and Monte Carlo. Results are visualised using matplotlib, by plotting position against intensity for 1d, and including an intensity colour gradient for 2d.

## Requirements
- Python 3.x
- numpy
- scipy
- matplotlib
- sys
- time

## How to use
1. Run `moon_probe.py`.
2. Choose a timescale (long ~1 lunar year or short ~3-4 days).
3. The script will plot Moon and probe trajectories
4. Automatically checks energy conservation.

## Example Output

Here is an example simulation of one lunar year for the Moon and probe:

![Moon and Probe Trajectories](lunar_year_figure.png)

## Features
- 2D Plots of Moon and probe orbits 
- Interactive menu for timescale 
- Energy conservation check

## Author
James Barnes – Physics Student, University of Bristol
