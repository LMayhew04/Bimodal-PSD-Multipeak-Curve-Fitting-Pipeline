README – Particle Size Modelling
================================

Overview:
---------
This MATLAB pipeline fits particle size distributions using models with
2 or 3 peaks. It automatically selects the best model using the BIC score.

NEW FEATURES:
-------------
✓ Automatically chooses between 2- and 3-component fits.
✓ Reports and displays residual statistics (Max Residual, RMSE).
✓ Annotates each peak on the plot with its diameter.
✓ Optional bootstrapping for confidence intervals.
✓ Supports weighted fitting using adjustable exponent (alpha).

Customizing Input Columns:
--------------------------
Default column names expected in Excel:
  - 'ParticleDiameter'
  - 'SizeDistribution1'

If your file uses different names:
Edit this in 'runAllFits.m':
    customNames = {'YourDiameterColumn', 'YourDistributionColumn'}

If your file has no headers:
Edit loadParticleData call like this:
    'DiameterIdx', 6, 'DistIdx', 8

Fitting Options:
----------------
In 'runAllFits.m' you can control:
  - sheetToUse = 'all'   % or 0 for all but 1st, or specific numbers
  - weightExponent = 0.5 % affects emphasis on small diameters
  - useBootstrap = true  % enables resampling for confidence intervals

Output:
-------
- Figures with fit + annotations saved in: /Figures/
- Fit summary CSV with model, BIC, and residual metrics: Fit_Summary.csv
- Optionally: confidence interval summaries via bootstrapping
