#!/usr/bin/env python3

"""
CABLE site run with full biogeochemistry (CNP) and POP
=======================================================

- Spin-up: using pre-industrial CO2, NDEP, PDEP. Currently this is using
           values from AmazonFACE experiment, i.e. CO2=284.7;
           NDEP-0.79 kg N ha-1 yr-1; PDEP=0.144 kg P ha-1 yr-1
- Transient: 1851-1998, varying CO2, NDEP and PDEP with recycled met data
- Simulation: actual experiment CO2, NDEP/PDEP & met

Options to turn use C/CN/CNP with POP on/off

That's all folks.

"""

__author__ = "Martin De Kauwe, Vanessa Haverd"
__version__ = "1.0 (21.07.2018)"
__email__ = "mdekauwe@gmail.com, Vanessa.Haverd@csiro.au"

import os
import sys
import glob
import shutil
import tempfile
import pandas as pd
import xarray as xr
import numpy as np

class RunCable(object):

    def __init__(self, experiment_id, namelist_dir, param_dir, output_dir, restart_dir,
                 dump_dir, met_fname, co2_ndep_fname, nml_fn, site_nml_fn,
                 veg_param_fn,log_dir, exe, aux_dir, biogeochem, call_pop,
                 verbose):

        self.experiment_id = experiment_id
        self.namelist_dir = namelist_dir
        self.param_dir = param_dir
        self.output_dir = output_dir
        self.restart_dir = restart_dir
        self.dump_dir = dump_dir
        self.met_fname = met_fname
        self.co2_ndep_fname = co2_ndep_fname
        self.nml_fn = nml_fn
        self.site_nml_fn = site_nml_fn
        self.veg_param_fn = veg_param_fn
        self.cable_restart_fname = "%s_cable_rst.nc" % (self.experiment_id)
        self.casa_restart_fname = "%s_casa_rst.nc" % (self.experiment_id)
        self.pop_restart_fname = "%s_pop_rst.nc" % (self.experiment_id)
        self.climate_restart_fname = "%s_climate_rst.nc" % (self.experiment_id)
        self.log_dir = log_dir
        self.cable_exe = exe
        self.aux_dir = aux_dir
        self.verbose = verbose
        self.nyear_spinup = 30
        self.limit_labile =".FALSE."
        if biogeochem == "C":
            self.biogeochem = 1
            #self.vcmax = "standard"
            self.vcmax = "Walker2014"
            self.vcmax_feedback = ".TRUE."
        elif biogeochem == "CN":
            self.biogeochem = 2
            self.vcmax = "Walker2014"
            self.vcmax_feedback = ".TRUE."
        elif biogeochem == "CNP":
            self.biogeochem = 3
            self.vcmax = "Walker2014"
            self.vcmax_feedback = ".TRUE."
        else:
            raise ValueError("Unknown biogeochemistry option: C, CN, CNP")
        self.call_pop = call_pop
        if self.call_pop:
            self.pop_flag = ".TRUE."
        else:
            self.pop_flag = ".FALSE."
        self.debug = True

    def main(self, SPIN_UP=False, TRANSIENT=False, SIMULATION=False):

        num = 1
        not_in_equilibrium = True

        (st_yr, en_yr,
         st_yr_trans, en_yr_trans,
         st_yr_spin, en_yr_spin) = self.get_years()

        self.initial_setup(st_yr_spin, en_yr_spin, st_yr, en_yr)

        if SPIN_UP == True:

            # initial spin
            print("Spinup")
            self.limit_labile =".TRUE."
            self.run_me()
            self.clean_up(end=False, tag="zero")
            #sys.exit()

            while num < 4:
                self.logfile="log_ccp%d" % (num)
                self.setup_re_spin(number=num)
                self.run_me()
                self.clean_up(end=False, tag="ccp%d" % (num))
                self.logfile="log_sa%d" % (num)
                self.setup_analytical_spin(number=num, st_yr_spin=st_yr_spin,
                                           en_yr_spin=en_yr_spin )
                self.run_me()
                self.clean_up(end=False, tag="saa%d" % (num))
                not_in_equilibrium = self.check_steady_state(num)
                num += 1
                
            not_in_equilibrium=True
            self.limit_labile =".FALSE."
            while not_in_equilibrium:

                self.logfile="log_ccp%d" % (num)
                self.setup_re_spin(number=num)
                self.run_me()
                self.clean_up(end=False, tag="ccp%d" % (num))
                #sys.exit()
                self.logfile="log_sa%d" % (num)
                self.setup_analytical_spin(number=num, st_yr_spin=st_yr_spin,
                                           en_yr_spin=en_yr_spin )
                self.run_me()
                self.clean_up(end=False, tag="saa%d" % (num))
                #sys.exit()
                not_in_equilibrium = self.check_steady_state(num)

                num += 1

            # one final spin
            self.logfile="log_ccp%d" % (num)
            self.setup_re_spin(number=num)
            self.run_me()
            self.clean_up(end=False, tag="ccp%d" % (num))

        if TRANSIENT == True:
            print("Transient")

            self.setup_transient(st_yr_trans, en_yr_trans, st_yr, en_yr)
            self.run_me()
            self.clean_up(end=False, tag="transient")

        if SIMULATION == True:
            print("Simulation")

            self.setup_simulation(st_yr, en_yr)
            self.run_me()

        self.clean_up(end=True)

    def get_years(self):
        """
        Figure out the start and end of the met file, the number of times we
        need to recycle the met data to cover the transient period and the
        start and end of the transient period.
        """
        pre_indust = 1850

        ds = xr.open_dataset(self.met_fname)

        st_yr = pd.to_datetime(ds.time[0].values).year

        # PALS met files final year tag only has a single 30 min, so need to
        # end at the previous year, which is the real file end
        en_yr = pd.to_datetime(ds.time[-1].values).year - 1

        # length of met record
        nrec = en_yr - st_yr + 1

        # number of times met data is recycled during transient simulation
        nloop_transient = np.ceil((st_yr - 1 - pre_indust) / nrec) - 1

        # number of times met data is recycled with a spinup run of nyear_spinup
        nloop_spin = np.ceil( self.nyear_spinup / nrec)

        st_yr_transient = st_yr - 1 - nloop_transient * nrec + 1
        en_yr_transient = st_yr_transient + nloop_transient * nrec - 1

        en_yr_spin = st_yr_transient - 1
        st_yr_spin = en_yr_spin - nloop_spin * nrec + 1

        return (st_yr, en_yr, st_yr_transient, en_yr_transient,
                st_yr_spin, en_yr_spin)

    def adjust_nml_file(self, fname, replacements):
        """
        Adjust the params/flags in the CABLE namelise file. Note this writes
        over whatever file it is given!

        Parameters:
        ----------
        fname : string
            parameter filename to be changed.
        replacements : dictionary
            dictionary of replacement values.

        """
        f = open(fname, 'r')
        param_str = f.read()
        f.close()
        new_str = self.replace_keys(param_str, replacements)
        fd, path = tempfile.mkstemp()
        os.write(fd, str.encode(new_str))
        os.close(fd)
        shutil.copy(path, fname)
        os.remove(path)

    def replace_keys(self, text, replacements_dict):
        """ Function expects to find CABLE namelist file formatted key = value.

        Parameters:
        ----------
        text : string
            input file data.
        replacements_dict : dictionary
            dictionary of replacement values.

        Returns:
        --------
        new_text : string
            input file with replacement values

        """
        lines = text.splitlines()
        for i, row in enumerate(lines):
            # skip blank lines
            if not row.strip():
                continue
            if "=" not in row:
                lines[i] = row
                continue
            elif not row.startswith("&"):
                key = row.split("=")[0]
                val = row.split("=")[1]
                lines[i] = " ".join((key.rstrip(), "=",
                                     replacements_dict.get(key.strip(),
                                     val.lstrip())))

        return '\n'.join(lines) + '\n'

    def initial_setup(self, st_yr_spin, en_yr_spin, st_yr, en_yr):
        """
        Setup CABLE namelist file for spinup from zero
        """
        shutil.copyfile(os.path.join(self.namelist_dir, "site.nml"),
                        self.site_nml_fn)
        shutil.copyfile(os.path.join(self.namelist_dir, "cable.nml"),
                        self.nml_fn)
        #self.add_missing_options_to_nml_file(self.nml_fn)

        out_fname = os.path.join(self.output_dir,
                                 "%s_out_cable_zero.nc" % (self.experiment_id))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        out_fname_CASA = os.path.join(self.output_dir,
                                 "%s_out_CASA_zero.nc" % (self.experiment_id))
        if os.path.isfile(out_fname_CASA):
            os.remove(out_fname_CASA)

        out_fname_POP = os.path.join(self.output_dir,
                                 "%s_out_POP_zero.nc" % (self.experiment_id))
        if os.path.isfile(out_fname_POP):
            os.remove(out_fname_POP)

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_zero.txt" % (self.experiment_id))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        replace_dict = {
                        "RunType": '"spinup"',
                        "CO2NDepFile": "'%s'" % (self.co2_ndep_fname),
                        "spinstartyear": "%d" % (st_yr),
                        "spinendyear": "%d" % (en_yr),
                        "spinCO2": "284.7",
                        "spinNdep": "0.79",
                        "spinPdep": "0.144",
        }
        self.adjust_nml_file(self.site_nml_fn, replace_dict)

        replace_dict = {
                        "filename%met": "'%s'" % (self.met_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "casafile%out": "'%s'" % (out_fname_CASA),
                        "cable_user%POP_outfile" : "'%s'" % (out_fname_POP),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_out": "'%s'" % os.path.join(self.restart_dir, self.cable_restart_fname),
                        "cable_user%climate_restart_out": "'%s'" % os.path.join(self.restart_dir, self.climate_restart_fname),
                        "cable_user%POP_restart_out": "'%s'" % os.path.join(self.restart_dir, self.pop_restart_fname),
                        "casafile%cnpepool": "'%s'" % os.path.join(self.restart_dir, self.casa_restart_fname),
                        "filename%restart_in": "''" ,
                        "cable_user%climate_restart_in": "''" ,
                        "cable_user%POP_restart_in": "''",
                        "filename%type": "'%s'" % (os.path.join(self.aux_dir, "offline/gridinfo_CSIRO_1x1.nc")),
                        "filename%veg": "'%s'" % os.path.join(self.param_dir, veg_param_fn),
                        "filename%soil": "'%s'" % os.path.join(self.param_dir, soil_param_fn),
                        "output%restart": ".TRUE.",
                        "casafile%phen": "'%s'" % (os.path.join(self.aux_dir, "core/biogeochem/modis_phenology_csiro.txt")),
                        "casafile%cnpbiome": "'%s'" % (os.path.join(self.param_dir, bgc_param_fn)),
                        "cable_user%RunIden": "'%s'" % (self.experiment_id),
                        "cable_user%POP_out": "'rst'",
                        "cable_user%POP_rst": "'./'",
                        "cable_user%POP_fromZero": ".T.",
                        "cable_user%CASA_fromZero": ".T.",
                        "cable_user%CLIMATE_fromZero": ".T.",
                        "cable_user%vcmax": "'%s'" % (self.vcmax),
                        "cable_user%YearStart": "%d" % (st_yr_spin),
                        "cable_user%YearEnd": "%d" % (en_yr_spin),
                        "cable_user%CASA_SPIN_STARTYEAR": "%d" % (st_yr_spin),
                        "cable_user%CASA_SPIN_ENDYEAR": "%d" % (en_yr_spin),
                        "cable_user%CALL_POP": "%s" % (self.pop_flag),
                        "output%averaging": "'monthly'",
                        "icycle": "%d" % (self.biogeochem),
                        "l_vcmaxFeedbk": "%s" % (self.vcmax_feedback),
                        "l_laiFeedbk": ".FALSE.",         
                        "cable_user%CASA_OUT_FREQ":  "'annually'" ,
                        "cable_user%limit_labile": "%s" % (self.limit_labile) ,
                        "cable_user%SRF":  ".T." ,  
                        "cable_user%STRF": "'DAMM'" ,
                        "cable_user%SMRF": "'DAMM'"
          }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_re_spin(self, number=None):
        """
        Adjust the CABLE namelist file with the various flags for another spin
        """
        out_log_fname = "%s_log_ccp%d.txt" % (self.experiment_id, number)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        out_fname = "%s_out_cable_ccp%d.nc" % (self.experiment_id, number)
        out_fname = os.path.join(self.output_dir, out_fname)
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        out_fname_CASA = "%s_out_CASA_ccp%d.nc" % (self.experiment_id, number)
        out_fname_CASA = os.path.join(self.output_dir, out_fname_CASA)
        if os.path.isfile(out_fname_CASA):
            os.remove(out_fname_CASA)

        out_fname_POP = "%s_out_POP_ccp%d.nc" % (self.experiment_id, number)
        out_fname_POP = os.path.join(self.output_dir, out_fname_POP)
        if os.path.isfile(out_fname_POP):
            os.remove(out_fname_POP)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_in": "'%s'" % os.path.join(self.restart_dir, self.cable_restart_fname),
                        "cable_user%climate_restart_in": "'%s'" % os.path.join(self.restart_dir, self.climate_restart_fname),
                        "cable_user%POP_restart_in": "'%s'" % os.path.join(self.restart_dir, self.pop_restart_fname),
                        "casafile%cnpipool": "'%s'" % os.path.join(self.restart_dir,self.casa_restart_fname),
                        "cable_user%POP_fromZero": ".F.",
                        "cable_user%CASA_fromZero": ".F.",
                        "cable_user%POP_rst": "'./'",
                        "cable_user%CLIMATE_fromZero": ".F.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".TRUE.",
                        "cable_user%CASA_NREP": "0",
                        "cable_user%SOIL_STRUC": "'sli'",
                        "icycle": "%d" % (self.biogeochem),
                        "leaps": ".TRUE.",
                        "spincasa": ".FALSE.",
                        "casafile%c2cdumppath": "' '",
                        "output%restart": ".TRUE.",
                        "filename%out": "'%s'" % (out_fname),
                        "casafile%out": "'%s'" % (out_fname_CASA),
                        "cable_user%POP_outfile" : "'%s'" % (out_fname_POP),
                        "cable_user%limit_labile": "%s" % (self.limit_labile) ,
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_analytical_spin(self, number, st_yr_spin, en_yr_spin):
        """
        Adjust the CABLE namelist file with the various flags for the
        analytical spin step
        """
        out_log_fname = "%s_log_analytic_%d.txt" % (self.experiment_id, number)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        out_fname_CASA = "%s_out_CASA_analytic_%d.nc" % \
                            (self.experiment_id, number)
        out_fname_CASA = os.path.join(self.output_dir, out_fname_CASA)
        if os.path.isfile(out_fname_CASA):
            os.remove(out_fname_CASA)

        out_fname_POP = "%s_out_POP_analytic_%d.nc" % (self.experiment_id, number)
        out_fname_POP = os.path.join(self.output_dir, out_fname_POP)
        if os.path.isfile(out_fname_POP):
            os.remove(out_fname_POP)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "icycle": "%d" % (self.biogeochem + 10), # Need to add 10 for spinup
                        "cable_user%CASA_DUMP_READ": ".TRUE.",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "cable_user%CASA_NREP": "1",
                        "cable_user%SOIL_STRUC": "'default'",
                        "leaps": ".FALSE.",
                        "spincasa": ".TRUE.",
                        "casafile%c2cdumppath": "'./'",
                        "cable_user%CASA_SPIN_STARTYEAR": "%d" % (st_yr_spin),
                        "cable_user%CASA_SPIN_ENDYEAR": "%d" % (en_yr_spin),
                        "casafile%out": "'%s'" % (out_fname_CASA),
                        "cable_user%POP_outfile" : "'%s'" % (out_fname_POP),
                        "cable_user%limit_labile": "%s" % (self.limit_labile) ,
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_transient(self, st_yr_trans, en_yr_trans, st_yr, en_yr):
        """
        Adjust the CABLE namelist file for the transient run, i.e. 1850 to XXXX
        """
        replace_dict = {
                        "RunType": '"transient"',
                        "CO2NDepFile": "'%s'" % (self.co2_ndep_fname),
                        "spinstartyear": "%d" % (st_yr),
                        "spinendyear": "%d" % (en_yr),
           }
        self.adjust_nml_file(self.site_nml_fn, replace_dict)

        out_log_fname = "%s_log_transient.txt" % (self.experiment_id)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        out_fname = "%s_out_cable_transient.nc" % (self.experiment_id)
        out_fname = os.path.join(self.output_dir, out_fname)
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        out_fname_CASA = "%s_out_casa_transient.nc" % (self.experiment_id)
        out_fname_CASA = os.path.join(self.output_dir, out_fname_CASA)
        if os.path.isfile(out_fname_CASA):
            os.remove(out_fname_CASA)

        replace_dict = {
                        "filename%out": "'%s'" % (out_fname),
                        "filename%log": "'%s'" % (out_log_fname),
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "cable_user%CASA_NREP": "0",
                        "cable_user%SOIL_STRUC": "'sli'",
                        "output%restart": ".TRUE.",
                        "output%averaging": "'monthly'",
                        "spinup": ".FALSE.",
                        "icycle": "%d" % (self.biogeochem),
                        "POPLUC": ".T.",
                        "filename%out": "'%s'" % (out_fname),
                        "casafile%out": "'%s'" % (out_fname_CASA),
                        "cable_user%YearStart": "%d" % (st_yr_trans),
                        "cable_user%YearEnd": "%d" % (en_yr_trans),
                        "filename%restart_in": "'%s'" % os.path.join(self.restart_dir, self.cable_restart_fname),
                        "cable_user%climate_restart_in": "'%s'" % os.path.join(self.restart_dir, self.climate_restart_fname),
                        "cable_user%POP_restart_in": "'%s'" % os.path.join(self.restart_dir, self.pop_restart_fname),
                        "casafile%cnpipool": "'%s'" % os.path.join(self.restart_dir,self.casa_restart_fname),
                        "cable_user%POP_fromZero": ".F.",
                        "cable_user%CASA_fromZero": ".F.",
                        "cable_user%CLIMATE_fromZero": ".F.",
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_simulation(self, st_yr, en_yr):
        """
        Adjust the CABLE namelist file for the experiment years
        """
        replace_dict = {
                        "RunType": '"historical"',
                        "CO2NDepFile": "'%s'" % (self.co2_ndep_fname),
                        "spinstartyear": "%d" % (st_yr),
                        "spinendyear": "%d" % (en_yr)
        }
        self.adjust_nml_file(self.site_nml_fn, replace_dict)

        out_log_fname = "%s_log_simulation.txt" % (self.experiment_id)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        out_fname = "%s_out_cable.nc" % (self.experiment_id)
        out_fname = os.path.join(self.output_dir, out_fname)

        out_fname_CASA = "%s_out_casa.nc" % (self.experiment_id)
        out_fname_CASA = os.path.join(self.output_dir, out_fname_CASA)

        if os.path.isfile(out_fname):
            os.remove(out_fname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "output%averaging": "'daily'",
                        "icycle": "%d" % (self.biogeochem),
                        "cable_user%YearStart": "%d" % (st_yr),
                        "cable_user%YearEnd": "%d" % (en_yr),
                        "filename%out": "'%s'" % (out_fname),
                        "POPLUC": ".F.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "cable_user%SOIL_STRUC": "'sli'",
                        "spincasa": ".FALSE.",
                        "spinup": ".FALSE.",
                        "output%averaging": "'all'",
                        "casafile%out": "'%s'" % (out_fname_CASA),
                        "filename%restart_in": "'%s'" % os.path.join(self.restart_dir, self.cable_restart_fname),
                        "cable_user%climate_restart_in": "'%s'" % os.path.join(self.restart_dir, self.climate_restart_fname),
                        "cable_user%POP_restart_in": "'%s'" % os.path.join(self.restart_dir, self.pop_restart_fname),
                        "casafile%cnpipool": "'%s'" % os.path.join(self.restart_dir,self.casa_restart_fname),
                        "cable_user%POP_fromZero": ".F.",
                        "cable_user%CASA_fromZero": ".F.",
                        "cable_user%CLIMATE_fromZero": ".F.",
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def run_me(self):
        if self.verbose:
            os.system("%s" % (self.cable_exe))
        else:
            # No outputs to the screen, stout and stderr to dev/null
            os.system("%s > /dev/null 2>&1" % (self.cable_exe))

    def check_steady_state(self, num):
        """
        Check whether the plant (leaves, wood and roots) and soil
        (fast, slow and active) carbon pools have reached equilibrium. To do
        this we are checking the state of the last year in the previous spin
        cycle to the state in the final year of the current spin cycle.
        """
        tol = 0.2 # This is quite high, I use 0.005 in GDAY
        g_2_kg = 0.001

        if num == 1:
            prev_cplant = 99999.9
            prev_csoil = 99999.9
        else:
            fname = "%s_out_CASA_ccp%d.nc" % (self.experiment_id, num-1)
            fname = os.path.join(self.output_dir, fname)
            ds = xr.open_dataset(fname)
            prev_cplant = ds.cplant[:,:,0].values[-1].sum() * g_2_kg
            prev_csoil = ds.csoil[:,:,0].values[-1].sum() * g_2_kg

        fname = "%s_out_CASA_ccp%d.nc" % (self.experiment_id, num)
        fname = os.path.join(self.output_dir, fname)
        ds = xr.open_dataset(fname)
        new_cplant = ds.cplant[:,:,0].values[-1].sum() * g_2_kg
        new_csoil = ds.csoil[:,:,0].values[-1].sum() * g_2_kg

        if ( np.fabs(prev_cplant - new_cplant) < tol and
             np.fabs(prev_csoil - new_csoil) < tol ):
            not_in_equilibrium = False
        else:
            not_in_equilibrium = True

        if self.debug:
            print("*", num, not_in_equilibrium,
                  "*cplant", np.fabs(prev_cplant - new_cplant),
                  "*csoil", np.fabs(prev_csoil - new_csoil))

        return not_in_equilibrium

    def clean_up(self, end=True, tag=None):
        """
        Move restart files to a directory and delete various files we no longer
        need that CABLE spits out as it spins up.
        """
        if end:
            for f in glob.glob("c2c_*_dump.nc"):
                shutil.move(f, os.path.join(self.dump_dir, f))
            f = "cnpfluxOut.csv"
            if os.path.isfile(f):
                os.remove(f)
            f = "new_sumbal"
            if os.path.isfile(f):
                os.remove(f)
            for f in glob.glob("*.out"):
                os.remove(f)
            for f in glob.glob("restart_*.nc"):
                os.remove(f)
        else:
            old = os.path.join(self.restart_dir, self.cable_restart_fname)
            new = "%s_%s.nc" % (old[:-3], tag)
            shutil.copyfile(old, new)

            old = os.path.join(self.restart_dir, self.casa_restart_fname)
            new = "%s_%s.nc" % (old[:-3], tag)
            shutil.copyfile(old, new)

            old = os.path.join(self.restart_dir, self.climate_restart_fname)
            new = "%s_%s.nc" % (old[:-3], tag)
            shutil.copyfile(old, new)

            if self.call_pop:
                old = os.path.join(self.restart_dir, self.pop_restart_fname)
                new = "%s_%s.nc" % (old[:-3], tag)
                shutil.copyfile(old, new)

    def add_missing_options_to_nml_file(self, fname, line_start=None):
        """
        Some of the flags we may wish to change are missing from the default
        file so we can't adjust them via this script...add them
        """
        if line_start is None:
            line_start = sum(1 for line in open(fname)) - 1

        f = open(fname, "r")
        contents = f.readlines()
        f.close()

        arg = "   cable_user%GW_MODEL = .FALSE.\n"
        contents.insert(line_start, arg)
        line_start += 1

        arg = "   cable_user%or_evap = .FALSE.\n"
        contents.insert(line_start, arg)
        line_start += 1

        arg = "   filename%gw_elev = 'GSWP3_elevation_slope_stddev.nc'\n"
        contents.insert(line_start, arg)
        line_start += 1

        tmp_fname = "tmp.nml"
        f = open(tmp_fname, "w")
        contents = "".join(contents)
        f.write(contents)
        f.close()

        shutil.move(tmp_fname, fname)


if __name__ == "__main__":

    cwd = os.getcwd()
    namelist_dir = "namelists"
    param_dir = "params"
    dump_dir = "dump"
    met_dir = "met"
    co2_ndep_dir = "met"
    aux_dir = "./CABLE-AUX/"
    log_dir = "logs"
    output_dir = "outputs"
    restart_dir = "restart_files"
    nml_fn = "cable.nml"
    site_nml_fn = "site.nml"
    #met_fname = os.path.join(met_dir, '%s.1.4_met.nc' % (experiment_id))
    met_fname = os.path.join(met_dir, 'AU_Cum_2014_2017_met_LAIv.nc')
    co2_ndep_fname = os.path.join(co2_ndep_dir,
                                  "AmaFACE_co2npdepforcing_1850_2100_AMB.csv")
    veg_param_fn = "Cumberland_veg_params.txt"
    bgc_param_fn = "Cumberland_pftlookup.csv"
    soil_param_fn = "def_soil_params.txt"   # only used when soilparmnew = .FALSE. in cable.nml
    exe = "./cable"

    # special for PEST
    #veg_param_fn = "def_veg_params_pest.txt"
    #bgc_param_fn = "pftlookup_pest.csv"
    #param_dir = "./"



    call_pop = True
    verbose = False

    if not os.path.exists(restart_dir):
        os.makedirs(restart_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    if not os.path.exists(dump_dir):
        os.makedirs(dump_dir)

    #for biogeochem in ["C", "CN", "CNP"]:
    for biogeochem in ["CNP"]:

        experiment_id = "Cumberland_POP_%s" % (biogeochem)
        C = RunCable(experiment_id, namelist_dir, param_dir, output_dir, restart_dir,
                     dump_dir, met_fname, co2_ndep_fname, nml_fn, site_nml_fn,
                     veg_param_fn, log_dir, exe, aux_dir, biogeochem, call_pop,
                     verbose)
        C.main(SPIN_UP=True, TRANSIENT=True, SIMULATION=True)
