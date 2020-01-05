# dispersioneering
A Matlab dispersion relation solver for magnetizied plasmas (both cold and hot).

# Octave notes
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
