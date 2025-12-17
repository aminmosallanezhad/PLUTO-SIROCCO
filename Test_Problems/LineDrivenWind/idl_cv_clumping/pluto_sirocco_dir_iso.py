#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Script: pluto_sirocco_dir_iso.py
# Author: Nick Higginbottom
# Revised By: Amin Mosallanezhad (a.mosallanezhad@soton.ac.uk)
# Last Modified: 2024-11-19
# -----------------------------------------------------------------------------

import subprocess
import sys
import glob
import os
import shutil
from astropy import constants as c
from astropy import units as u
from scipy.integrate import quad
import pyPLUTO.pload as pp
import numpy as np
import pluto_sirocco_sub as pss  # Importing custom module for PLUTO-Python subroutines
import re
from pluto_sirocco_config import data  # Import configuration data

# -------------------------------------------------------------------------
# Logging setup
# -------------------------------------------------------------------------

# Open a log file to record the run progress as early as possible.
# buffering=1 → line buffered (flushes at newlines); we also flush manually.
logfile = open("run_logfile", "w", buffering=1)


def log(message: str) -> None:
    """
    Log a message to both stdout and run_logfile, always with a newline and flush.
    """
    if not message.endswith("\n"):
        out = message + "\n"
    else:
        out = message
    # Write to logfile
    logfile.write(out)
    logfile.flush()
    # Echo to stdout
    sys.stdout.write(out)
    sys.stdout.flush()


# -------------------------------------------------------------------------
# Code units from PLUTO
# -------------------------------------------------------------------------

# Attempt to import code units from definitions.h using the get_units function from pluto_sirocco_sub.
# This is necessary because PLUTO expects geometrical parameters in code units,
# and the translation between code and physical units is set by the UNIT_ commands.
try:
    UNIT_DENSITY, UNIT_LENGTH, UNIT_VELOCITY = pss.get_units()
    log("Successfully read code units from definitions.h")
except Exception as e:
    log(f"Unable to open definitions.h file - big problems ({e})")
    logfile.close()
    sys.exit(1)

# -------------------------------------------------------------------------
# Initial parameters for the simulation run
# -------------------------------------------------------------------------

t0 = 1.256637061  # Initial run time for the PLUTO simulation (in seconds)
dt = 1.256637061  # Time increment between calls to PLUTO (in seconds)

# Ensure t0 is at least 1.0
if t0 == 0.0:
    log("We need to run for at least one second, dummy")
    t0 = 1.256637061

init_sir_cycles = 2   # Initial number of cycles that SIROCCO will perform
sir_cycles = 2        # Number of cycles in SIROCCO once started

# -------------------------------------------------------------------------
# Handle restarts by checking for command-line arguments
# -------------------------------------------------------------------------

if len(sys.argv) < 2:
    # No restart argument provided; start from the beginning
    istart = 0
    sim_time = t0
    log("No restart cycle specified; starting from cycle 0")
else:
    # Attempt to parse the restart cycle number from the command-line argument
    try:
        istart = int(sys.argv[1])
    except ValueError:
        log(f"Invalid argument: {sys.argv[1]} is not an integer.")
        logfile.close()
        sys.exit(1)

    if istart > 0:
        # Restarting from a previous cycle
        root = f"{istart:08d}"
        directory = "cycle" + root
        log(f"We will be trying to restart using files in {directory}")

        # Copy necessary files from the previous cycle's directory
        try:
            for filename in os.listdir(directory):
                source_file = os.path.join(directory, filename)
                if os.path.isfile(source_file):
                    shutil.copy(source_file, ".")
        except Exception as e:
            log(f"Cannot restart: {e}")
            logfile.close()
            sys.exit(1)

        # Update simulation time based on the last PLUTO output
        try:
            last_pluto = pp.pload(istart)
            log(f"Last run finished at {last_pluto.SimTime}")
            sim_time = last_pluto.SimTime + dt
        except Exception as e:
            log(f"Error loading PLUTO restart data for cycle {istart}: {e}")
            logfile.close()
            sys.exit(1)
    else:
        # Starting from the first cycle
        istart = 0
        sim_time = t0
        log("Restart index <= 0; starting from cycle 0")

# -------------------------------------------------------------------------
# In the absence of force multiplier files, use a k-alpha formulation
# -------------------------------------------------------------------------

data["k"] = 0.59
data["alpha"] = -0.6  # Alpha is normally negative

# Copy initial flux files needed for the simulation
try:
    shutil.copy("directional_flux_theta_init.dat", "directional_flux_theta.dat")
    shutil.copy("directional_flux_phi_init.dat", "directional_flux_phi.dat")
    shutil.copy("directional_flux_r_init.dat", "directional_flux_r.dat")
    log("Copied initial directional flux files")
except Exception as e:
    log(f"Error copying directional flux files: {e}")
    logfile.close()
    sys.exit(1)

# -------------------------------------------------------------------------
# Main simulation loop
# -------------------------------------------------------------------------

for i in range(istart, 1000):
    # For cycles beyond the initial run, set k and alpha to dummy values
    if i > 0:
        data["k"] = 999
        data["alpha"] = 999

    # Format cycle number for consistent file naming
    root = f"{i:08d}"
    log(f"Making a PLUTO input file for cycle {i}")

    # Create a PLUTO input file using the current simulation time and data
    pss.pluto_input_file(sim_time, data)

    # Determine the command to run PLUTO, considering the number of processors
    if i == 0:
        # First cycle; no restart
        if data["nproc_pluto"] == 1:
            cmd = ["./pluto"]
        else:
            cmd = ["mpirun", "-np", str(data["nproc_pluto"]), "./pluto"]
    else:
        # Subsequent cycles; restart from the previous cycle
        if data["nproc_pluto"] == 1:
            cmd = ["./pluto", "-restart", str(i)]
        else:
            cmd = ["mpirun", "-np", str(data["nproc_pluto"]), "./pluto", "-restart", str(i)]

    log("Running PLUTO run")
    log("Command line: " + " ".join(cmd))

    # Execute the PLUTO simulation
    try:
        with open("pluto_log", "w") as log_file:
            result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            log("Error: PLUTO command failed. Check pluto_log for details.")
            logfile.close()
            sys.exit(1)
    except Exception as e:
        log(f"Failed to run PLUTO: {e}")
        logfile.close()
        sys.exit(1)

    log("Finished PLUTO run")

    # Update simulation time and cycle number for the next iteration
    ifile = i + 1
    sim_time = sim_time + dt
    root = f"{ifile:08d}"
    dbl = f"data.{ifile:04d}.dbl"
    directory = "cycle" + root

    log(f"Turning dbl file {dbl} into a SIROCCO model file")

    # Convert PLUTO output to a SIROCCO-readable model file
    try:
        sir_model_file = pss.pluto2sir_rtheta(ifile)
    except Exception as e:
        log(f"Error converting PLUTO output to SIROCCO model file for cycle {ifile}: {e}")
        logfile.close()
        sys.exit(1)

    log(f"Made a SIROCCO model file called {sir_model_file}")

    # Create a SIROCCO parameter file for the simulation
    log("Making a SIROCCO input file")
    pss.sirocco_input_file(root, data, cycles=init_sir_cycles + i * sir_cycles)
    log("Successfully made a SIROCCO input file")

    # Copy the SIROCCO parameter file to a generic name (input.pf) to maintain consistency
    try:
        shutil.copy(f"{root}.pf", "input.pf")
    except Exception as e:
        log(f"Error copying {root}.pf to input.pf: {e}")
        logfile.close()
        sys.exit(1)

    log(f"Command line: cp {root}.pf input.pf")

    # Running the SIROCCO simulation
    if i > 0:
        # For cycles beyond the initial run, restart SIROCCO using the previous windsave file
        log("Restarting SIROCCO")

        # Update densities in the windsave file using the new model
        log("Running modify_wind on the old windsave")
        cmd = ["modify_wind" + data["RAD_CODE_VER"], "-model_file", sir_model_file, "input"]
        log("Command line: " + " ".join(cmd))

        # Execute modify_wind to update the windsave file
        try:
            with open("mod_wind_log", "w") as log_file:
                result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
            if result.returncode != 0:
                log("Error: modify_wind command failed. Check mod_wind_log for details.")
                logfile.close()
                sys.exit(1)
        except Exception as e:
            log(f"Failed to run modify_wind: {e}")
            logfile.close()
            sys.exit(1)

        # Copy the new windsave file to input.wind_save
        try:
            shutil.copy("new.wind_save", "input.wind_save")
        except Exception as e:
            log(f"Error copying new.wind_save to input.wind_save: {e}")
            logfile.close()
            sys.exit(1)

        log("Copied new.wind_save to input.wind_save")

        # Prepare the command to run SIROCCO in restart mode
        cmd = [
            "mpirun",
            "-np",
            str(data["nproc_sir"]),
            data["RAD_CODE"] + data["RAD_CODE_VER"],
            "-f",
            "-r",
            "input.pf",
        ]
    else:
        # First run of SIROCCO; no restart
        cmd = [
            "mpirun",
            "-np",
            str(data["nproc_sir"]),
            data["RAD_CODE"] + data["RAD_CODE_VER"],
            "-f",
            "input.pf",
        ]

    log("Running SIROCCO")
    log("Command line: " + " ".join(cmd))

    # Execute the SIROCCO simulation
    try:
        with open("sirocco_log", "w") as log_file:
            result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            log("Error: SIROCCO command failed. Check sirocco_log for details.")
            logfile.close()
            sys.exit(1)
    except Exception as e:
        log(f"Failed to run SIROCCO: {e}")
        logfile.close()
        sys.exit(1)

    log("Finished SIROCCO")

    # ---------------------------------------------------------------------
    # Process SIROCCO outputs to generate files needed by PLUTO
    # ---------------------------------------------------------------------

    cmd = [
        "mpirun",
        "-np",
        str(data["nproc_radhydro"]),
        "rad_hydro_files" + data["RAD_CODE_VER"],
        "input",
    ]

    log("Running rad_hydro_files")
    log("Command line: " + " ".join(cmd))

    try:
        with open("rad_hydro_files_output", "w") as log_file:
            # run, but DO NOT stop if return code is non-zero
            result = subprocess.run(
                cmd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
            )

        # Instead of exiting on error → only warn
        if result.returncode != 0:
            warning_msg = (
                f"WARNING: rad_hydro_files returned exit code {result.returncode}.\n"
                "Ignoring error and continuing. See rad_hydro_files_output for details."
            )
            log(warning_msg)
    except Exception as e:
        # Also ignore unexpected exceptions, but log them
        warning_msg = (
            f"WARNING: rad_hydro_files failed to launch: {e}\n"
            "Continuing anyway."
        )
        log(warning_msg)

    # ---------------------------------------------------------------------
    # Copy the resulting Python data files with the current root name
    # ---------------------------------------------------------------------

    try:
        for src, dst in [
            ("py_heatcool.dat", f"{root}_py_heatcool.dat"),
            ("py_driving.dat", f"{root}_py_driving.dat"),
            ("py_pcon_data.dat", f"{root}_py_pcon_data.dat"),
            ("py_ion_data.dat", f"{root}_py_ion_data.dat"),
            ("py_spec_data.dat", f"{root}_py_spec_data.dat"),
            ("py_fluxes.dat", f"{root}_py_fluxes.dat"),
        ]:
            if os.path.exists(dst):
                os.remove(dst)
            shutil.copy(src, dst)
    except Exception as e:
        log(f"Error copying py_* files: {e}")
        logfile.close()
        sys.exit(1)

    # ---------------------------------------------------------------------
    # Run the external CAK (Castor-Abbott-Klein) code
    # ---------------------------------------------------------------------

    cmd = ["mpirun", "-np", str(data["nproc_cak"]), "./cak_v3"]
    log("Running CAK")
    log("Command line: " + " ".join(cmd))

    try:
        with open("cak_output", "w") as log_file:
            result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            log("Error: cak_v3 command failed. Check cak_output for details.")
            logfile.close()
            sys.exit(1)
    except Exception as e:
        log(f"Failed to run cak_v3: {e}")
        logfile.close()
        sys.exit(1)

    log("Finished CAK")

    # ---------------------------------------------------------------------
    # Create prefactors file for the Blondin heating and cooling rates
    # ---------------------------------------------------------------------

    log("Creating prefactors file for the Blondin heating and cooling rates")
    try:
        pss.pre_calc(ifile)
    except Exception as e:
        log(f"Error creating prefactors file in pre_calc for cycle {ifile}: {e}")
        logfile.close()
        sys.exit(1)
    log("Finished creating prefactors file")

    # ---------------------------------------------------------------------
    # Clean up old directories to manage storage space
    # ---------------------------------------------------------------------

    old_directory = f"cycle{ifile - 10:08d}"
    if os.path.exists(old_directory):
        try:
            shutil.rmtree(old_directory)
            log(f"Removed old directory {old_directory}")
        except Exception as e:
            log(f"Error removing old directory {old_directory}: {e}")

    # Create a new directory for the current cycle to store output files
    try:
        os.mkdir(directory)
    except FileExistsError:
        try:
            shutil.rmtree(directory + "_old")
        except FileNotFoundError:
            pass
        shutil.move(directory, directory + "_old")
        os.mkdir(directory)

    # Copy all relevant files into the cycle directory for storage
    try:
        for filename in [
            "dbl.out",
            "pluto.ini",
            "restart.out",
            "grid.out",
            "definitions.h",
            dbl,
        ]:
            dst = os.path.join(directory, os.path.basename(filename))
            if os.path.exists(dst):
                os.remove(dst)
            shutil.copy(filename, directory)

        # Copy all Python data files matching the pattern py_*.dat
        for file in glob.glob("py_*.dat"):
            dst = os.path.join(directory, os.path.basename(file))
            if os.path.exists(dst):
                os.remove(dst)
            shutil.copy(file, directory)

        for filename in [
            "M_UV_data.dat",
            "prefactors.dat",
            "directional_flux_r.dat",
            "directional_flux_theta.dat",
            "directional_flux_phi.dat",
            "input.wind_save",
        ]:
            dst = os.path.join(directory, os.path.basename(filename))
            if os.path.exists(dst):
                os.remove(dst)
            shutil.copy(filename, directory)

        # Move the model and parameter files into the cycle directory
        for filename in [
            "py_pcon_data.dat",
            sir_model_file,
            f"{root}.pf",
            f"{root}_py_heatcool.dat",
            f"{root}_py_driving.dat",
            f"{root}_py_pcon_data.dat",
            f"{root}_py_ion_data.dat",
            f"{root}_py_spec_data.dat",
            f"{root}_py_fluxes.dat",
        ]:
            dst = os.path.join(directory, os.path.basename(filename))
            if os.path.exists(dst):
                os.remove(dst)
            shutil.move(filename, directory)

    except Exception as e:
        log(f"Error copying files to {directory}: {e}")
        logfile.close()
        sys.exit(1)

    log("Finished tidying up")

# -------------------------------------------------------------------------
# End of run
# -------------------------------------------------------------------------

logfile.close()
print("Fin")
