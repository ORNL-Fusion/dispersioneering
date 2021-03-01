# dispersioneering
A Matlab dispersion relation solver for magnetizied plasmas (both cold and hot). It will always solve the cold dispersion via the quadratic method, but will then addionally add solutions from a root solver with either the hot or cold dielectrics. 

# Quickstart
```
cd dispersioneering/
addpath(genpath('./'))
dispersioneering(@case1,'num_points',30,'use_root_finder',true,'use_cold_eps',false)
```

# Running a new case
The `dispersioneering/cases` folder contains the input files to specify the case you want to run. The way to specifiy the case you want is in the call to `dispersioneering()` in the first arguement, i.e., 

the `dispersioneering/cases/fast_wave.m` case is run using
```
dispersioneering(@fast_wave)
```

The case input files specify function handles for each of the backgroun profiles. 

## Options
`numpoints` is an integer specifying how many, equally spaced points to sample your profiles with. 
`use_root_finder` is a `true` or `false` (default) for using a root finder to solve for the determinant == 0 or not. If not it defaults to cold and uses the quadratic expression. This is applicable to hot or cold. 
`use_cold_eps` is a `true` or `false` (default) for selecting the hot or cold dielectric. 
`kper_min` is a numeric value specifying the minimum value of `kper` (default=-1000) to scan over for the root finder.
`kper_max` is a numeric value specifying the maximum value of `kper` (default=+1000) to scan over for the root finder.

# Limitations
- At present only solves for `kper` (given `kpar`), not `kpar` (given `kper`);
- At present assumes Maxwellian for hot plasma. 

# Dependencies
You will need the optimization matlab toolbox for hot plasma calculations. 


# Octave notes - just don't do it. 
Octave is extremely slow (~100x slower), and really only supported to allow for reviewers to run the code if they do not have access to Matlab with the optimization toolbox. 

## Dependencies

within Octave install the following (this can take a while)
```
pkg install -forge struct io statistics optim netcdf
```
then you can load the `optim` package
```
pkg load optim netcdf
```
