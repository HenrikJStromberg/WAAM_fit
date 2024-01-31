# WAAM_fit
Evaluates the fitness of stuctures for WAAM manufacturing.
The conceptual idea is to derive geometry indicators (radii, radii gradients and angles) by inscibing heuvers' speheres in the geometry.
Technically this is solved by applying a medial axis transform (MAT). Afterwards, the geometry indicators are mapped to geometry features and compared against design constraints.
The results are displayed using Gmsh where all non-manufacturable elements get highlighted.

## How to use
The main part of the application is a python program which is invoked from the command line with, while `<input-file>` is a geometry file in the `.stp`-format:
```
waam_fit <input-file> -o <output-file> -s <mesh-size>
```
See the examples dir for example input files.

## Installation
However, to run the python program, compilation and installation of the rust bindings and the command-line-tool is needed, which is handeld by maturin
It is recommended to install the tool in a python virtual envrionment (either `venv` or `conda`). 
In the virtual environment install the python dependencies specified in `requirements.txt`. Fo example with using `pip` call :
```
pip install -r requirements.txt
```
in the project directory.

Afterwards, the commands
```
maturin build --release
```
and
```
pip install waam_fit
``` 
have to be executed to compile the rust bindings and install the python script as command-line-tool (see also the [maturin package documentation](https://www.maturin.rs)).

## Configuration
The Post-Processing can be configured via the `WAAM.toml` file (in [/src/waam_fit/WAAM.toml](/src/waam_fit/WAAM.toml)). 

The `feature` group defines the mapping of geometry indicators to corresponding features. For every feature manufacturability is determined by the geometry indicator specified as `data`, 
while the `min`, `max` and `scale` attributes configure the data range that is not manufacturable. To only include elements belonging to the feature a `filter` can be specified.

A `filter` consists of another geometry indicator specified as `data` and a threshold defined by `less_eq` or `greater_eq` to configure the included elements (all elements with `data`$\leq$ / $\geq$ the threshold).
