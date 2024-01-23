# Uber Quick Guide

Assuming you have computing resources on gadi@NCI, installing and running CABLE is as simple as:

1. Clone the source code:

        git clone https://github.com/CABLE-LSM/CABLE.git
        cd CABLE

1. Build a serial version of CABLE

        ./build.bash

1. Execute this serial version of CABLE

        ./bin/cable

## In slightly more detail

The CABLE executable (./cable) is configured via the namelist cable.nml. By default the provided cable.nml points to 4-year meteorological forcing from the [Tumbarumba][fluxnet-tumba] fluxnet site. NB: spinning up the model has been switched off in cable.nml (spinup=.FALSE.)

We say the "serial" version as CABLE does support a parallel configuration as well.

ACCESS-ESM1.5 uses a version of CABLE based on CABLE-2.3.4.

ACCESS-CM2 uses a version of CABLE that is closer to the HEAD of the trunk at the time of writing (December 2019).

For a more detailed discussion of CABLE offline, please see the other sections of the [Cable's User Guide][userguide].

[fluxnet-tumba]: https://fluxnet.org/sites/siteinfo/AU-Tum
[userguide]: index.md
