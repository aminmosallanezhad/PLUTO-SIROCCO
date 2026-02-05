#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Script: pluto_sirocco_init.py
# Purpose: Initialize and run the SIROCCO simulation pipeline with PLUTO.
# Author: Amin Mosallanezhad, a.mosallanezhad@soton.ac.uk
# Date: 2025-12-17
#
# Key fix:
#   - rad_hydro_files now uses data["nproc_radhydro"] (fallback: nproc_sir, then 1)
#
# Notes:
#   - Defensive logging for missing optional files and external command failures
#   - Variables are declared first, assigned later (per requested code style)
# -----------------------------------------------------------------------------

import os
import sys
import shutil
import subprocess
from datetime import datetime

from pluto_sirocco_config import data
import pluto_sirocco_sub as pss


# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------

logfile_path = "processing_log.txt"
logfile = None


def log_and_print(message: str) -> None:
    """Log a message to both the terminal and logfile (line-buffered)."""
    global logfile
    print(message)
    sys.stdout.flush()
    if logfile is not None:
        logfile.write(message + "\n")
        logfile.flush()


def safe_remove(path: str) -> None:
    """Remove a file if it exists, ignore if it doesn't."""
    if os.path.exists(path):
        try:
            os.remove(path)
        except Exception as e:
            log_and_print(f"WARNING: Could not remove {path}: {e}")


def cleanup() -> None:
    """Clean up temporary files in case of errors."""
    files_to_remove = [
        "input.pf",
        "sirocco_log",
        "rad_hydro_files_output",
        "pluto_log",
    ]
    for f in files_to_remove:
        safe_remove(f)
    log_and_print("Cleanup finished.")


def run_command(
    cmd,
    stdout_path=None,
    check=False,
    ignore_fail=False,
    label="command",
):
    """
    Run a subprocess command with optional stdout redirection.
    - check=True raises CalledProcessError on non-zero exit.
    - ignore_fail=True logs warning and continues on non-zero exit.
    """
    result = None
    log_and_print(f"Executing {label}: {' '.join(cmd)}")

    try:
        if stdout_path is not None:
            with open(stdout_path, "w") as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    check=check,
                )
        else:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=check,
                text=True,
            )
            if result.stdout:
                snippet = result.stdout[-8000:]  # last 8k chars
                log_and_print(f"{label} output (tail):\n{snippet}")

        if result is not None and result.returncode != 0:
            msg = f"WARNING: {label} returned exit code {result.returncode}."
            if ignore_fail:
                log_and_print(msg + " Continuing anyway.")
            else:
                log_and_print(msg)

        return result

    except subprocess.CalledProcessError as e:
        log_and_print(f"ERROR: {label} failed with exit code {e.returncode}.")
        if ignore_fail:
            log_and_print(f"WARNING: Ignoring failure in {label} and continuing.")
            return None
        raise

    except FileNotFoundError as e:
        log_and_print(f"ERROR: Could not run {label} (missing executable?): {e}")
        if ignore_fail:
            log_and_print(f"WARNING: Ignoring launch failure in {label} and continuing.")
            return None
        raise

    except Exception as e:
        log_and_print(f"ERROR: Unexpected failure running {label}: {e}")
        if ignore_fail:
            log_and_print(f"WARNING: Ignoring unexpected failure in {label} and continuing.")
            return None
        raise


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> int:
    # --- declarations at top ---
    init_sir_cycles = None
    sim_time = None
    ifile = None
    root = None
    dbl = None

    source_pf = None
    dest_pf = None

    nproc_sir = None
    sir_exe = None
    sir_cmd = None

    nproc_radhydro = None
    radhydro_exe = None
    radhydro_cmd = None

    run_stamp = None

    # --- assignments below ---
    data["k"] = 0.59
    data["alpha"] = -0.6  # Alpha is normally negative

    init_sir_cycles = 30
    sim_time = 0.0
    ifile = 0
    root = f"{ifile:08d}"
    dbl = f"data.{ifile:04d}.dbl"

    run_stamp = datetime.now().strftime("%A, %Y-%m-%d %H:%M:%S")

    log_and_print("Starting SIROCCO processing pipeline")
    log_and_print("Script Run Date: " + run_stamp)

    # Resource monitoring
    nproc_sir = int(data.get("nproc_sir", 1))
    log_and_print(f"Using {nproc_sir} processor(s) for SIROCCO (nproc_sir).")

    # Step 1: Generate PLUTO input file
    log_and_print("Step 1: Generating PLUTO input file.")
    try:
        pss.pluto_input_file(sim_time, data)
    except Exception as e:
        log_and_print(f"ERROR: pluto_input_file failed: {e}")
        cleanup()
        return 1

    # Step 2: Run PLUTO executable (best effort)
    log_and_print("Step 2: Running PLUTO executable.")
    try:
        run_command(["./pluto"], stdout_path="pluto_log", check=False, ignore_fail=True, label="PLUTO")
    except Exception:
        pass

    # Step 3: Generate Python model file
    log_and_print(f"Step 3: Converting {dbl} into a Python model file.")
    try:
        sirocco_model_file = pss.pluto2sir_rtheta(ifile)
        log_and_print(f"Generated Python model file: {sirocco_model_file}")
    except Exception as e:
        log_and_print(f"ERROR: pluto2sir_rtheta failed: {e}")
        cleanup()
        return 1

    # Step 4: Create SIROCCO parameter file (.pf)
    log_and_print("Step 4: Creating SIROCCO input parameter file (.pf).")
    try:
        pss.sirocco_input_file(root, data, cycles=init_sir_cycles)
        log_and_print("SIROCCO parameter file created successfully.")
    except Exception as e:
        log_and_print(f"ERROR: sirocco_input_file failed: {e}")
        cleanup()
        return 1

    # Step 5: Copy parameter file to input.pf
    log_and_print("Step 5: Copying parameter file to input.pf.")
    source_pf = f"{root}.pf"
    dest_pf = "input.pf"

    if not os.path.exists(source_pf):
        log_and_print(f"ERROR: Source file {source_pf} does not exist.")
        cleanup()
        return 1

    try:
        shutil.copy(source_pf, dest_pf)
        log_and_print(f"Successfully copied {source_pf} to {dest_pf}")
    except Exception as e:
        log_and_print(f"ERROR copying {source_pf} to {dest_pf}: {e}")
        cleanup()
        return 1

    # Step 6: Run SIROCCO with mpirun (strict: fail if SIROCCO fails)
    log_and_print("Step 6: Running SIROCCO simulation with mpirun.")
    try:
        sir_base = data.get("RAD_CODE", "sirocco")
        sir_ver = data.get("RAD_CODE_VER", "")
        sir_exe = f"{sir_base}{sir_ver}"

        sir_cmd = [
            "mpirun",
            "-np",
            str(nproc_sir),
            sir_exe,
            "-f",
            "-p",
            "2",
            "input.pf",
        ]
        run_command(sir_cmd, stdout_path="sirocco_log", check=True, ignore_fail=False, label="SIROCCO")
    except Exception as e:
        log_and_print(f"ERROR: SIROCCO failed: {e}")
        cleanup()
        return 1

    # Step 7: Run rad_hydro_files (FIXED mpirun -np)
    log_and_print("Step 7: Running rad_hydro_files and saving output.")

    try:
        nproc_radhydro = int(data.get("nproc_radhydro", nproc_sir if nproc_sir else 1))
    except Exception:
        nproc_radhydro = nproc_sir if nproc_sir else 1

    log_and_print(f"Using {nproc_radhydro} processor(s) for rad_hydro_files (nproc_radhydro).")

    radhydro_exe = "rad_hydro_files" + data.get("RAD_CODE_VER", "")
    radhydro_cmd = [
        "mpirun",
        "-np",
        str(nproc_radhydro),
        radhydro_exe,
        "input",
    ]

    # Warn & continue if rad_hydro_files fails (matching your full pipeline behavior)
    try:
        run_command(
            radhydro_cmd,
            stdout_path="rad_hydro_files_output",
            check=False,
            ignore_fail=True,
            label="rad_hydro_files",
        )
    except Exception:
        pass

    # Step 8: Rename flux files (best effort)
    log_and_print("Step 8: Renaming flux files.")
    rename_pairs = [
        ("directional_flux_r.dat", "directional_flux_r_init.dat"),
        ("directional_flux_theta.dat", "directional_flux_theta_init.dat"),
        ("directional_flux_phi.dat", "directional_flux_phi_init.dat"),
    ]

    for src, dst in rename_pairs:
        if not os.path.exists(src):
            log_and_print(f"WARNING: {src} not found; skipping rename to {dst}.")
            continue
        try:
            if os.path.exists(dst):
                safe_remove(dst)
            shutil.move(src, dst)
            log_and_print(f"Renamed {src} -> {dst}")
        except Exception as e:
            log_and_print(f"WARNING: Could not rename {src} -> {dst}: {e}")

    log_and_print("Completed SIROCCO run and file management successfully.")
    return 0


if __name__ == "__main__":
    # Open logfile as early as possible (line-buffered)
    try:
        logfile = open(logfile_path, "w", buffering=1)
    except Exception:
        logfile = None

    try:
        rc = main()
    finally:
        if logfile is not None:
            logfile.close()

    sys.exit(rc)
