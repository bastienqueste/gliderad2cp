---
title: 'gliderad2cp: A Python package to process Nortek AD2CP velocity profiles from gliders'

tags:
  - Python
  - oceanography
  - currents
  - adcp
  - shear

authors:
  - name: Bastien Y. Queste
    orcid: 0000-0002-3786-2275
    corresponding: true
    equal-contrib: true
    affiliation: "1, 2"
  - name: Callum Rollo
    orcid: 0000-0002-5134-7886
    equal-contrib: true
    affiliation: 2
  - name: Estel Font
    orcid: 0000-0002-8767-6665
    affiliation: 1

affiliations:
 - name: Department of Marine Science, University of Gothenburg, Natrium, Box 463, 405 30 Göteborg, Sweden
   index: 1
 - name: Voice of the Ocean Foundation, Skeppet Ärans väg 3, 426 71 Västra Frölunda, Sweden 
   index: 2

date: 18 October 2023

bibliography: paper.bib
---

# Summary

Oceanographers routinely measure ocean currents to understand and map the transport of ocean properties. Measuring currents is most commonly done using instrument called acoustic doppler current profilers (ADCP). These instruments emit shorts pings of sound and listen for the echoing soundwaves which bounce off of water molecules and suspended particles. These return echoes contain much valuable information. The delay between emission and receiving tells us distance to the particles, and the pitch change of the echo tells us the relative velocity of the particles to the sensor. Using multiple beams of sounds, the ADCP can determine 3-dimensional currents at range. ADCP are however limited by power and size; there is a direct trade-off between size, power and transducer capability. There is a also a trade-off between ping frequency and effective range before the soudn wave is attenuated. Large ocean going vessels can carry low frequency ADCP with ranges of hundreds of meters down into the water column, but with resolutions in the tens of meters.

Ocean gliders are small, low power, autonomous underwater vehicles which glider up and down in the water column, collecting measurements of ocean properties throughout. Ocean gliders now have the ability to carry small ADCP such as the Nortek Glider AD2CP, with 4 beams and a frequency of 1MHz. The high frequency means that the sensor can only measure currents up to approximately 15 m away from the glider; however as the glider travels up and down through the water column, coverage is possible down to the glider's full depth. The key difficulty comes as the ADCP measure ocean currents relative to the glider, rather than relative to ground. As the glider's velocity is an order of magnitude greater than ocean current velocities in most areas, a different form of processing is required known as the lowered ADCP method. This toolbox collects successive measurements of ocean currents as the glider profiles up and down and performs the following steps, proividing figures for easy assessment of processing quality:
- Clean the ADCP data and remove bad measurements.
- Correct the vertical alignment (in the earth frame of reference) of velocity measurements across all beams.
- Convert the velocity data from ADCP-relative (*ie*. beam direction; Fig. \autoref{fig:beam2xyz}), to glider-relative (*ie*. X, Y, Z) and finally to earth-relative velocities (*ie*. East, North, Up).
- Calculate the vertical gradient in earth-relative velocities, also known as vertical shear (*ie*. $$\frac{\delta v}{\delta z}$$).
- Reconstruct full-depth profiles of vertical shear from the successive low-range measurements to small scale relative changes in ocean currents, but lacking an absolute reference.
- Determine the mean ocean current over the period of the glider dive by comparing ADCP-derived glider speed through water to its GPS-derived speed over land, the difference being caused by ocean currents.
- Reference the full high-resolution vertical shear profile using the glider's dive-averaged current to provide a high-resolution absolute measurements of ocean currents.

![ADCP beams measure the along-beam velocity which needs to be converted to X,Y,Z velocities relative to teh glider's frame of reference. The coordinate transform matrix is specific to each instrument as it is defined by the angle of the different beams relative to the glider.\label{fig:beam2xyz}](paper_figures/beam2xyz.png)

# Statement of need

Software for processing ADCP data exists, with tools provided by instrument manufacturers, private companies and open-source communities. However, none of these tools apply the lowered ADCP method on Nortek Glider AD2CP sensors. 

Similar for Teledyne (https://github.com/JGradone/Glider_ADCP_Real_Time_Processing), (https://ieeexplore.ieee.org/document/7098134)

A great portion of the processing can be automated in a straightforward way and require minimal user input; namely the successive rounds of 3-dimensional coordinate transforms and regridding of velocity data along isobars which may deviate due to glider pitch. This toolbox greatly simplifies the file handling, integration of glider data to ADCP data, and complex trigonometry necessary to obtain high quality shear data.

Furthermore, the integration of the Nortek AD2CP varies across glider manufacturers, either using alternating 3-beam configurations between up and down profiles (on the Seaglider or the Spray) or using 4 beams at all times (on the SeaExplorer). This python package allows users to easily load Nortek AD2CP netCDF files and pull the raw data to provide clean shear estimates with consistent processing and quality control independent of which glider they use.

This toolbox integrates work performed at GU, VOTO, but also Tanaka and Todd while basing on best practices.

@Callum : need to integrate new OceanGlider format as input to make it nicely platform independent. What's update on new format specification?

# Package description

## Quality control

- Velocity

- Amplitude

- Correlation

- Issues with first bin and side lobe interference

## Coordinate transformations

- Regridding to avoid shear smearing.

Shear smearing

![The Nortek AD2CP measurements are time-gated at the same intervals for each individual beam, meaning that the relation between echo delay and measurement range is the same for all 4 beams and does not account for the more open front and back beam angles. The purpose is to have 3 beams at equal angles from vertical when the glider is diving at the correct angle (17.4$$^\circ$$ from horizontal for the Nortek AD2CP; in grey on the left). If the glider is flying at a different angle, there will be a mismatch in depth between the 3 beams (in gray on the right) which requires regridding and use of different bins (in green on the right) to minimise shear smearing.\label{fig:regridding}](paper_figures/regridding.png)

- Standard matrices for the Nortek AD2CP.

@Callum: should we add functionality to extract coordinate transform matrix from netCDF metadata or keep it hard coded?

## Integration with glider data

- What is necessary format of glider data, why do we require these variables

- Lat and lon data needed where/when and need for accurate GPS data?

# Usage

## Assessing shear data quality

- Example code

- Figure list and description

## Obtaining referenced velocities


![L-ADCP method.\label{fig:ladcp}](paper_figures/lADCP.png)

- Built in DVL approach to calculate DAC and reference

- Using glider flight model

- Using fixed point reference

- Using bottom tracking

- Visbeck LSQ approach to multiple constraints

# Known issues

- Shear bias

- Compass calibrations

Figures can be included like this:

![Caption for example figure.\label{fig:example}](paper_figures/blank.png)

and referenced from text using \autoref{fig:example}.

# Final words

# Citations

# Acknowledgements

ONR-Global grant

Formas grant

Voice of the Ocean Foundation

Intern students from the French Ecole Navale (for figures)

# References
