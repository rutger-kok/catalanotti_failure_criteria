# 3D Failure Criteria for Composite Materials
Implementation of the 3D failure criteria for fiber reinforced composites developed by Catalanotti et al. (2013) in fixed format Fortran 77.

## Summary
This repository contains implements the three-dimensional failure criteria for fiber reinforced composites developed by Catalanotti et al. [1]. Required in-situ material properties can be calculated using the insitu_prop_calculator.py script. Failure envelopes can be plotted using the create_envelopes.py script. To call the Fortran scripts from Python you will need to use F2PY.

Alternatively, the subroutines can be used as a part of an Abaqus user subroutine, or called from a standalone Fortran program.

## Usage
To use the subroutines as a part of an Abaqus simulation your installation must be linked with a Fortran compiler and compatible Visual Studio installation, see: 

The following material properties must be defined to plot the failure envelopes:   
**X<sub>T</sub>** = tensile strength fiber direction (GPa)
**X<sub>C</sub>** = compressive strength fiber direction (GPa)
**Y<sub>T</sub><sup>is</sup>** =  in-situ tensile strength transverse direction (GPa)
**Y<sub>C</sub><sup>is</sup>** = in-situ compressive strength transverse direction (GPa)
**S<sub>L</sub><sup>is</sup>** = in-situ longitudinal shear strength (GPa)
**η<sub>L</sub>** = shear friction coefficient longitudinal direction (-)
**α<sub>0</sub>** = failure plane angle pure transverse compression (degrees)

If the in-situ properties are unknown, they can be calculated using the insitu_prop_calculator.py script. This will require the following material properties
**β** = shear response factor
**G<sub>12</sub>** = shear modulus (MPa)
**t** = cured ply thickness (mm)
**E<sub>11</sub>** = longitudinal modulus (MPa)
**E<sub>22</sub>** = transverse modulus (MPa)
**ν<sub>21</sub>** =  Poisson's ratio 21-direction (-)
**G<sub>Ic</sub>** = mode I fracture toughness (N/mm)
**G<sub>IIc</sub>** = mode II fracture toughness (N/mm)

## List of Fortran source code
- **catalanotti.f** : Catalanotti failure criteria

## List of Python source code
- **insitu_prop_calculator.py** : In-situ property calculator
- **failure_envelopes.py** : Plots failure envelopes according to Catalanotti failure criteria

***
Rutger Kok  
PhD Candidate  
email : rutger.kok@ed.ac.uk  

Institute for Infrastructure and Environment  
University of Edinburgh    
Thomas Bayes Road, King's Buildings, Edinburgh EH9 3FG   
United Kingdom

***
>[1] G. Catalanotti, P.P. Camanho, A.T. Marques  
>Three-dimensional failure criteria for fiber-reinforced laminates   
>Composite Structures 95 (2013) 63–79  
>http://dx.doi.org/10.1016/j.compstruct.2012.07.016  
