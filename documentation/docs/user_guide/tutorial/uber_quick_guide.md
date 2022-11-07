= Uber Quick Guide =

Assuming you have computing resources on gadi@NCI, installing and running CABLE is as simple as:

1. Checkout the source code:
{{{
svn checkout https://trac.nci.org.au/svn/cable/trunk MyCABLE
cd MyCABLE/offline
}}}

2. Build a serial version of CABLE
{{{
./build3.sh
}}}

3. Execute this serial version of CABLE
{{{
./cable
}}}


=== In slightly more detail ===

The CABLE executable (./cable) is configured via the namelist cable.nml. By default the provided cable.nml points to 4-year meteorological forcing from the [http://sites.fluxdata.org/AU-Tum/ Tumbarumba] fluxnet site. NB: spinning up the model has been switched off in cable.nml (spinup=.FALSE.)

We say the "serial" version as CABLE does support a  parralel configuration as well.

ACCESS-ESM1.5 uses  a version of CABLE based on CABLE-2.3.4.

ACCESS-CM2  uses  a version of CABLE that is closer to the HEAD of the trunk at the time of writing (December 2019).

!!For a more detailed discussion of CABLE offline please see [https://trac.nci.org.au/trac/cable/wiki/CableUserGuide this page] from the User Guide.
