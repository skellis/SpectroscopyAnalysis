iterateFitResFFT.py is the main operating file.
You can optionally load the real data included in this file (RGTCNQTMB.txt is the Raman Gain Matrix as a function of time delay 'timepointsTCNQTMB' in fs and Raman shift 'shiftxTCNQTMB' in cm-1 or if you would rather nm 'nmshiftxTCNQTMB') or you can generate dummyData which simulates a two dimensional excited state femtosecond stimulated Raman spectrum to operate on.
See https://pubs.acs.org/doi/abs/10.1021/jp5041986

The operating file calls methods in init.py, fitfuncs.py and graphing.py

init.py does the hard work of the iterative fitting. exctracting and sorting the residuals ect.

graphing.py has some plotting functions that are optially called at the bottom of iterateFitResFFT.py and should help trouble shoot future fittings.

the heavy numerical work is done by curve_fit which is finiky if you try to extend its functionality. I ran into some difficulty when I tried to get curve_fit to hold parameters constant. holdCurveFit.py and utilities.py work together to demonstrate a simple example of how this can be done by generating the correct curve_fit fit function as a string and then evaluating it with eval. I think this functionality will be invaluable for future scientists and statiscians using python.
