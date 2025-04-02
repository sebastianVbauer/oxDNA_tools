# oxDNA_tools

This is a collection of small Python scripts that can be used to analyze data generated from oxDNA simulations.

## Purposes
Quick overview:
* contour_length.py -- analysis tool to obtain contour lengths from oxDNA trajectories for arbitrarily complicated paths.
* forcex.py -- analysis tool to get force-extension curves from oxDNA trajectories and force files with two harmonic traps.
* forcys.py -- just like forcex.py, but for complete hysteresis (cf. reverse_traps.py).
* reverse_traps.py -- creates force files with two harmonic traps after force-extension simulations to reverse them.

## Usage
They are ready to be executed from the command line, taking required data as arguments. No adjustments in the code should be necessary.
Instructions on how to use the tools is given as introductory comments at the top of each file.

## Licence
Published under cc-by-sa 4.0, cf. https://creativecommons.org/licenses/by-sa/4.0/

-- Sebastian V. Bauer, sebastian.bauer@uni-mainz.de

-- Walther Lab at the University of Mainz, https://www.walther-group.com/
