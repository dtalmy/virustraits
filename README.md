justfe.py - MAIN code for fitting host-virus population dynamics to time-series data

ehuxbasic.f - fortran code that is used to generate a python wrapper with 'f2py -c ehuxbasic.f -m ehuxbasic' from the command line. This creates the module ehuxbasic that can be imported in any python script, and the function within it called

basicforshow.f - this just allows you to rerun the code at many time points (not just the sampling time points).. this is better for pplotting with
