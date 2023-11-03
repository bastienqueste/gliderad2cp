---
title: 'gliderad2cp Dcoumentation'

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

# Notes and TODO

Suggested reviewers: Laur Ferris, ?

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