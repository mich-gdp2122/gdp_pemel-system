# PEM Electrolyser Thermal Management System model

FEEG6013 Group Design Project, year 2021-2022

**Group 19, *Thermal Management of a Green Hydrogen Production Unit***

Created by Michael and Victor in MATLAB/Simulink


## Key Files & Folders:

### Files:

* `GDPPEMELSystem.prj`: Main project file; open this first.

* `MAIN.m`: Main program file; runs Simulink simulations and generates output plots and tables.

* `parameters_straight.m`, `parameters_serpentine.m` Place electrolyser input values in here, for straight and serpentine cooling channel geometries, respectively.

* `parameters.m`: Place other input parameters in here.

* `parameters_derived.m`: Calculates other required inputs, derived from above parameter files. No input parameters are placed in here.

* `HXsizer_*.m`: Heat-exchanger sizing scripts; calculates minimum required surface area for specified hot- and cold-side inlet/outlet temperatures.

* `ORCspec.m`: Calculates Organic Rankine Cycle (ORC) performance as well as thermodynamic properties at each cycle state.

* `data_*.m`: Calculates thermodynamic properties of corresponding fluid at given input conditions, using look-up tables from NIST Chemistry WebBook ([source](https://webbook.nist.gov/chemistry/)) and built-in Simscape fluid property components.

* `const_file.m`: Physical constants are placed in here.

* `output_*.mat`: Simulation output files automatically generated after running simulations; used by `MAIN.m` for generating output plots and tables.


### Folders:

* `TMS/`: Contains Simulink models of each systems configuration.

* `calculations/`: Contains calculation procedures and other functions used by multiple .m files.

* `labelpoints/`: MATLAB community add-on used in some plots ([source](https://uk.mathworks.com/matlabcentral/fileexchange/46891-labelpoints)).

* `resources/project`: Misc. project files.



## Directions:

1. Open project by opening `GDPPEMELSystem.prj` file. If unable to (eg. .prj file extension is used by a different application), open it instead by running the `_OPEN_PROJECT_Win10.bat` file (if using Windows) or `_OPEN_PROJECT_MacOS.sh` file (if using MacOS or Linux).

2. If necessary, check and/or modify any input values in the relevant parameter files described above (open them within the project).

3. Open and run the `MAIN.m` file within the project to run the model simulations and generate output results (will be saved into an automatically created output folder).
