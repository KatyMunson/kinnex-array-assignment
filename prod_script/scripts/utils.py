#!/usr/bin/env python3
"""
utils.py

Shared utilities for the Kinnex Iso-seq production pipeline.

Covers:
  - Manifest parsing (4- or 5-column format)
  - Arrays file parsing and KN_BC subsetting
  - Processing mode resolution (single-library vs multi-library),
    including auto-detection and manifest override with cross-validation
  - Assignment file loading (skips provenance comment headers)

Processing modes
----------------
  single-library : One library per KN_BC. Skera output is passed directly
                   through as HIGH_CONF. No Lima pass 1, no Bayesian
                   assignment, no split step.

  multi-library  : Two or more libraries share a KN_BC. Lima pass 1 feeds
                   the Bayesian classifier (assign_kinnex.py), which splits
                   ZMWs into HIGH_CONF / LOW_CONF / UNASSIGNED before Lima
                   pass 2 and IsoSeq refine.

Mode resolution order
---------------------
  1. If the manifest has a 5th column (mode override), use it — after
     cross-validating against the arrays file (see resolve_mode()).
  2. Otherwise, auto-detect from the number of libraries that share the
     sample's KN_BC in the arrays file. The decision is always logged.

Author:  KM
Created: 2026-02
"""

import os
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SINGLE = "single"   # single-library processing mode
MULTI  = "multi"    # multi-library (array assignment) processing mode

# ---------------------------------------------------------------------------
# Arrays file
# ---------------------------------------------------------------------------

def read_arrays(arrays_file):
    """
    Parse an arrays.txt file and return a dict mapping library name to its
    metadata.

    Each line in the file is tab-separated with no header:
        library_name  KN_BC  barcode1  barcode2  ...

    Lines with fewer than 3 fields are silently skipped (e.g. blank lines).

    Returns
    -------
    dict
        { library_name: {"kinnex": KN_BC, "barcodes": set(barcode, ...)} }
    """
    arrays = {}
    if not os.path.exists(arrays_file):
        raise FileNotFoundError(f"Arrays file not found: {arrays_file}")
    with open(arrays_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            lib_name = parts[0]
            kinnex   = parts[1]
            barcodes = set(parts[2:])
            arrays[lib_name] = {"kinnex": kinnex, "barcodes": barcodes}
    return arrays


def libraries_for_knbc(arrays, knbc):
    """
    Return a list of library names in `arrays` whose KN_BC matches `knbc`.

    Parameters
    ----------
    arrays : dict
        Output of read_arrays().
    knbc : str
        Kinnex barcode identifier (e.g. "bcM0003").

    Returns
    -------
    list of str
    """
    return [lib for lib, info in arrays.items() if info["kinnex"] == knbc]


# ---------------------------------------------------------------------------
# Mode resolution
# ---------------------------------------------------------------------------

def resolve_mode(sample, knbc, arrays_file, declared_mode=None):
    """
    Determine whether a sample should be processed in single-library or
    multi-library mode, and return the resolved mode string plus the list
    of libraries that will be processed.

    Auto-detection
    --------------
    If declared_mode is None, the number of libraries matching `knbc` in
    the arrays file determines the mode:
      - 1 library  → SINGLE
      - >1 library → MULTI
    The decision is printed to stdout so it is visible in the Snakemake log.

    Override cross-validation
    -------------------------
    If declared_mode is provided, it is honoured after cross-checking:

      declared MULTI, 1 library found
        → Warn and proceed as MULTI. The classifier will run against a single
          candidate; all ZMWs with enough segments will be called HIGH_CONF.

      declared SINGLE, >1 libraries found
        → Error. A library name cannot be chosen unambiguously from multiple
          candidates, so the pipeline cannot proceed.

      declared SINGLE, 1 library found
        → Proceed as SINGLE (explicit agreement with auto-detect).

      declared MULTI, >1 libraries found
        → Proceed as MULTI (explicit agreement with auto-detect).

    Parameters
    ----------
    sample : str
        Sample identifier (for logging / error messages).
    knbc : str
        Kinnex barcode for this sample (e.g. "bcM0003").
    arrays_file : str or Path
        Path to the arrays.txt file.
    declared_mode : str or None
        "single", "multi", or None (auto-detect).

    Returns
    -------
    tuple (mode, libraries)
        mode      : SINGLE or MULTI constant string
        libraries : list of library names that will be processed
    """
    arrays    = read_arrays(arrays_file)
    libs      = libraries_for_knbc(arrays, knbc)
    n_libs    = len(libs)

    if n_libs == 0:
        raise ValueError(
            f"Sample '{sample}': no libraries found for KN_BC '{knbc}' "
            f"in arrays file '{arrays_file}'. "
            f"Check that the KN_BC column matches exactly."
        )

    # ── No override: auto-detect ────────────────────────────────────────────
    if declared_mode is None:
        if n_libs == 1:
            mode = SINGLE
            print(
                f"[mode] {sample}: auto-detected SINGLE-LIBRARY "
                f"(1 library for KN_BC {knbc}: {libs[0]})",
                flush=True,
            )
        else:
            mode = MULTI
            print(
                f"[mode] {sample}: auto-detected MULTI-LIBRARY "
                f"({n_libs} libraries for KN_BC {knbc}: {', '.join(sorted(libs))})",
                flush=True,
            )
        return mode, libs

    # ── Override provided: cross-validate ──────────────────────────────────
    declared_mode = declared_mode.strip().lower()
    if declared_mode not in (SINGLE, MULTI):
        raise ValueError(
            f"Sample '{sample}': invalid mode '{declared_mode}' in manifest. "
            f"Must be 'single' or 'multi'."
        )

    if declared_mode == SINGLE and n_libs > 1:
        # Unresolvable: cannot pick a library name from multiple candidates.
        raise ValueError(
            f"Sample '{sample}': mode=single was declared in the manifest "
            f"but {n_libs} libraries were found for KN_BC '{knbc}' in "
            f"'{arrays_file}':\n"
            f"  {', '.join(sorted(libs))}\n"
            f"Cannot determine library name unambiguously. Either correct "
            f"the arrays file so only one library maps to '{knbc}', or "
            f"change the manifest mode to 'multi'."
        )

    if declared_mode == MULTI and n_libs == 1:
        # Unusual but valid: classifier will run against a single candidate.
        print(
            f"[mode] WARNING — {sample}: mode=multi was declared but only "
            f"1 library was found for KN_BC '{knbc}' ({libs[0]}). "
            f"Proceeding through the multi-library (array assignment) path "
            f"as declared. All ZMWs with sufficient segments will be called "
            f"HIGH_CONF.",
            flush=True,
        )
        return MULTI, libs

    # Declared mode agrees with the arrays file — proceed without noise.
    print(
        f"[mode] {sample}: mode={declared_mode.upper()}-LIBRARY "
        f"(declared; {n_libs} librar{'y' if n_libs == 1 else 'ies'} for "
        f"KN_BC {knbc})",
        flush=True,
    )
    return declared_mode, libs


# ---------------------------------------------------------------------------
# Manifest parsing
# ---------------------------------------------------------------------------

def parse_manifest(manifest_path, arrays_root=None):
    """
    Parse the production manifest and resolve processing modes.

    Manifest format (tab-separated, no required header, # lines ignored):

        sample   KN_BC   bam_path   arrays_path   [mode]

    The 5th column (mode) is optional. When absent, mode is auto-detected
    from the arrays file.

    Parameters
    ----------
    manifest_path : str or Path
        Path to the manifest file.
    arrays_root : str or Path or None
        If provided, relative arrays_path values in the manifest are resolved
        relative to this directory. Absolute paths are used as-is.

    Returns
    -------
    dict with keys:
        samples     : list of sample names in manifest order
        bams        : { sample: bam_path }
        array_files : { sample: arrays_path }
        knbcs       : { sample: KN_BC }
        modes       : { sample: SINGLE or MULTI }
        libraries   : { sample: [library_name, ...] }

    Raises
    ------
    FileNotFoundError
        If the manifest file does not exist.
    ValueError
        If any manifest line is malformed, or if mode cross-validation fails.
    """
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")

    samples     = []
    bams        = {}
    array_files = {}
    knbcs       = {}
    modes       = {}
    libraries   = {}

    with open(manifest_path) as f:
        for lineno, line in enumerate(f, 1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"Manifest line {lineno} has {len(fields)} field(s), "
                    f"expected at least 4 "
                    f"(sample, KN_BC, bam_path, arrays_path). "
                    f"Check for missing tabs or wrong delimiter.\n"
                    f"  Line content: {line.rstrip()!r}"
                )

            sample     = fields[0]
            knbc       = fields[1]
            bam        = fields[2]
            arrays_raw = fields[3]
            declared   = fields[4].strip() if len(fields) >= 5 else None

            # Resolve arrays path relative to arrays_root if needed
            if arrays_root and not os.path.isabs(arrays_raw):
                arrays_path = os.path.join(arrays_root, arrays_raw)
            else:
                arrays_path = arrays_raw

            if sample in samples:
                raise ValueError(
                    f"Manifest line {lineno}: duplicate sample name '{sample}'. "
                    f"Each sample must appear only once."
                )

            mode, libs = resolve_mode(sample, knbc, arrays_path, declared)

            samples.append(sample)
            bams[sample]        = bam
            array_files[sample] = arrays_path
            knbcs[sample]       = knbc
            modes[sample]       = mode
            libraries[sample]   = libs

    if not samples:
        raise ValueError(f"No samples found in manifest: {manifest_path}")

    return {
        "samples":     samples,
        "bams":        bams,
        "array_files": array_files,
        "knbcs":       knbcs,
        "modes":       modes,
        "libraries":   libraries,
    }


# ---------------------------------------------------------------------------
# Assignment file loading
# ---------------------------------------------------------------------------

def load_assignments(assign_path):
    """
    Load an assign_kinnex.py output file, skipping provenance comment lines
    and the column header.

    The assignment file format is tab-separated with a column header line
    and optional leading comment lines (starting with '#') that record
    provenance, run date, and scoring parameters:

        # assign_kinnex.py | 2026-02-14 12:00:00 | git: a3f2c1b
        # SPECIFIC_WEIGHT=1.0  MAX_SHARED_WEIGHT=0.2  ...
        ZMW  Assigned_Array  Classification  Top_Posterior  ...
        m64/1  LibA  HIGH_CONF  0.923  ...

    Parameters
    ----------
    assign_path : str or Path
        Path to the assignment .txt file.

    Returns
    -------
    list of dict
        One dict per ZMW row, keyed by column header names.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file contains no data rows after the header.
    """
    if not os.path.exists(assign_path):
        raise FileNotFoundError(f"Assignment file not found: {assign_path}")

    rows    = []
    header  = None

    with open(assign_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue          # skip blank lines and provenance comments
            if header is None:
                header = line.split("\t")
                continue          # first non-comment line is the column header
            fields = line.split("\t")
            rows.append(dict(zip(header, fields)))

    if header is not None and not rows:
        raise ValueError(
            f"Assignment file has a header but no data rows: {assign_path}"
        )

    return rows


def load_assignments_df(assign_path):
    """
    Like load_assignments(), but returns a pandas DataFrame.

    Requires pandas. Provenance comment lines (starting with '#') and the
    column header are handled automatically via pandas comment= parameter.

    Parameters
    ----------
    assign_path : str or Path

    Returns
    -------
    pandas.DataFrame
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError(
            "pandas is required for load_assignments_df(). "
            "Use load_assignments() for a dict-based alternative."
        )
    return pd.read_csv(assign_path, sep="\t", comment="#")


# ---------------------------------------------------------------------------
# Provenance header parsing
# ---------------------------------------------------------------------------

def parse_assignment_header(assign_path):
    """
    Parse the provenance comment lines from an assignment file header and
    return a dict of key-value pairs.

    Extracts run date, git hash, and all scoring parameters written by
    assign_kinnex.py. Useful for training scripts to warn when mixing
    assignment files generated with different parameters.

    Example header lines parsed:
        # assign_kinnex.py | 2026-02-14 12:00:00 | git: a3f2c1b-dirty
        # SPECIFIC_WEIGHT=1.0  MAX_SHARED_WEIGHT=0.2  DISCORDANT_PENALTY=-0.10
        # POSTERIOR_HIGH_CONF=0.840  POSTERIOR_LOW_CONF=0.50
        # MIN_OBS_HIGH_CONF=3  MIN_OBS_LOW_CONF=2

    Returns
    -------
    dict
        Keys include 'run_date', 'git', and any KEY=VALUE pairs found.
        Returns empty dict if no comment headers are present.
    """
    result = {}
    if not os.path.exists(assign_path):
        return result

    with open(assign_path) as f:
        for line in f:
            if not line.startswith("#"):
                break   # stop at first non-comment line (the header)
            line = line.lstrip("#").strip()
            if not line:
                continue

            # First header line: "assign_kinnex.py | DATE | git: HASH"
            if "assign_kinnex.py" in line:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) >= 2:
                    result["run_date"] = parts[1]
                if len(parts) >= 3 and parts[2].startswith("git:"):
                    result["git"] = parts[2].replace("git:", "").strip()
                continue

            # Subsequent lines: KEY=VALUE pairs separated by whitespace
            for token in line.split():
                if "=" in token:
                    k, v = token.split("=", 1)
                    result[k.strip()] = v.strip()

    return result
