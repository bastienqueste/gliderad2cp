gliderad2cp: Processing ad2cp data from gliders
=================================================

What is gliderad2cp?
---------------------

``gliderad2cp`` gliderad2cp processes data from the Nortek AD2CP acoustic doppler current profiler (ADCP) mounted on a glider. gliderad2cp takes data from the ADCP unit and the glider and combines them to produce estimates of vertical shear of velocity. It also prodives functionality to integrate these velocity shear profiles into absolute earth relative water vlocities.



Package description
-------------------------------

Quality control
^^^^^^^^^^^^^^^^^^^^^^
- Velocity

- Amplitude

- Correlation

- Issues with first bin and side lobe interference

Coordinate transformations
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Regridding to avoid shear smearing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Shear smearing

![The Nortek AD2CP measurements are time-gated at the same intervals for each individual beam, meaning that the relation between echo delay and measurement range is the same for all 4 beams and does not account for the more open front and back beam angles. The purpose is to have 3 beams at equal angles from vertical when the glider is diving at the correct angle (17.4$$^\circ$$ from horizontal for the Nortek AD2CP; in grey on the left). If the glider is flying at a different angle, there will be a mismatch in depth between the 3 beams (in gray on the right) which requires regridding and use of different bins (in green on the right) to minimise shear smearing.\label{fig:regridding}](paper_figures/regridding.png)

Standard matrices for the Nortek AD2CP.


Integration with glider data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- What is necessary format of glider data, why do we require these variables

- Lat and lon data needed where/when and need for accurate GPS data?

Usage
=========
Assessing shear data quality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Example code

- Figure list and description

Obtaining referenced velocities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

LADCP method
.. image::  ../../paper_figures/lADCP.png
![L-ADCP method.\label{fig:ladcp}](paper_figures/lADCP.png)

- Built in DVL approach to calculate DAC and reference

- Using glider flight model

- Using fixed point reference

- Using bottom tracking

- Visbeck LSQ approach to multiple constraints

Known issues
===============
- Shear bias

- Compass calibrations



.. toctree::
   :maxdepth: 3
   :caption: Contents:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
