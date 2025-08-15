#!/usr/bin/env python3
"""
Compute P_cb(k,z) with CAMB and (optionally) also total matter P(k) and transfer functions.

IMPORTANT CHANGE: Background tables (H(z) and optional Omega_cb(z)) are now written on a FIXED scale-factor grid
independent of the user-chosen P(k) redshift sampling. This avoids coupling background resolution to the (often
coarser) spectrum output list.

Grids:
  P(k) redshifts: controlled by --zmin/--zmax/--nz (must be >= 0 for spectra; negative z ignored for P(k)).
  Background grid: 200 linearly spaced scale factors a in [0.0001, 1.51356]. For each a we output z = 1/a - 1.
                   This spans very high redshift (z ~ 9999) down to z ~ -0.339 (a > 1 future times).

Outputs per P(k) redshift index i (always for z>=0):
    pk_*_cb_{i:03d}.dat      -> columns: k  P_cb  (prefix via --pk-prefix)

Additional outputs when --write-optional-files is set (for z>=0):
    pk_*_tot_{i:03d}.dat     -> total matter (cb + massive nu) spectrum
    tf_*_{i:03d}.dat         -> columns: k  T_cdm  T_baryon  T_photon  T_neutrino  T_nu  T_tot

Background outputs (independent grid):
    redshifts.dat            -> index  z_pk (only the P(k) list; unchanged)
    <hubble-file>            -> columns: z_bg  H(z_bg)   (fixed background grid)
    omega_cb.dat (optional)  -> columns: z_bg  Omega_cb(z_bg) (same background grid)

Unit control:
    --units h         : k in h/Mpc, P in (Mpc/h)^3 (default)
    --units physical  : k in 1/Mpc,  P in Mpc^3

Notes:
  - Power spectrum uses CAMB's interpolator var1=var2='delta_nonu' (CDM+baryons) and optionally 'delta_tot'.
  - Transfer functions: results.get_matter_transfer_data(); native k is k/h. Converted if physical units requested.
  - Transfer is log-interpolated onto the chosen k_grid.
  - CAMB is asked to compute matter power at the union of requested P(k) redshifts and the maximum background redshift
    so that H(z) is reliable up to z~9999.
"""

import argparse
import os
import numpy as np

import camb
from camb import model

def parse_redshifts(args) -> np.ndarray:
    """Return P(k) redshift list from user controls (may include negatives, but only z>=0 used for P(k))."""
    amin = 1.0 / (1 + args.zmax)
    amax = 1.0 / (1 + args.zmin)
    atab = np.linspace(amin, amax, args.nz, dtype=float)
    return 1.0 / atab - 1.0

def background_scale_factor_grid() -> tuple[np.ndarray, np.ndarray]:
    """Return (a_bg, z_bg) for the fixed background grid.

    a: 200 linearly spaced points in [1e-5, 1.5].
    z: 1/a - 1 for each of those scale factors.
    """
    a_bg = np.geomspace(1e-5, 1.5, 500, dtype=float)
    z_bg = 1.0 / a_bg - 1.0
    return a_bg, z_bg

def main():
    p = argparse.ArgumentParser(description="Compute CDM+baryons P(k) and transfer T(k) with CAMB.")
    # Redshift controls
    p.add_argument("--zmin", type=float, default=0.00, help="Minimum redshift (used with --zmax/--nz).")
    p.add_argument("--zmax", type=float, default=99.0, help="Maximum redshift (used with --zmin/--nz).")
    p.add_argument("--nz", type=int, default=100, help="Number of redshifts (used with --zmin/--zmax).")
    # Units choice
    p.add_argument("--units", choices=["h", "physical"], default="h",
                   help="Output units: 'h' => k[h/Mpc], P[(Mpc/h)^3]; 'physical' => k[1/Mpc], P[Mpc^3]. Default: h")
    # k-grid limits (interpreted according to --units)
    p.add_argument("--kmin", type=float, default=1e-4,
                   help="Minimum k in h/Mpc if --units h, or 1/Mpc if --units physical (default: 1e-4).")
    p.add_argument("--kmax", type=float, default=10.0,
                   help="Maximum k in h/Mpc if --units h, or 1/Mpc if --units physical (default: 10).")
    p.add_argument("--nk", type=int, default=512, help="Number of k points (log-spaced; default: 512).")
    # Cosmology (Planck-like defaults)
    p.add_argument("--H0", type=float, default=67.66, help="H0 in km/s/Mpc (default: 67.66).")
    p.add_argument("--Ob0", type=float, default=0.0489, help="Omega_b (default: 0.0489).")
    p.add_argument("--Ocdm0", type=float, default=0.2621, help="Omega_cdm (default: 0.2621).")
    p.add_argument("--ns", type=float, default=0.9649, help="Scalar spectral index n_s (default: 0.9649).")
    p.add_argument("--As", type=float, default=2.1e-9, help="Primordial amplitude A_s (default: 2.1e-9).")
    p.add_argument("--mnu", type=float, default=0.06, help="Sum of neutrino masses in eV (default: 0.06).")
    p.add_argument("--Neff", type=float, default=3.046, help="Effective N_eff for massless species (default: 3.046).")
    p.add_argument("--tau", type=float, default=0.054, help="Optical depth tau (default: 0.054).")
    # Linear / non-linear choice for P(k)
    p.add_argument("--nonlinear", action="store_true",
                   help="If set, include non-linear corrections for matter P(k). Linear by default.")
    # Output controls
    p.add_argument("--outdir", type=str, default=".", help="Output directory (default: current).")
    p.add_argument("--pk-prefix", type=str, default="pk", help="P(k) filename prefix (default: 'pk').")
    p.add_argument("--tf-prefix", type=str, default="tf", help="T(k) filename prefix (default: 'tf').")
    p.add_argument("--write-optional-files", action="store_true",
                   help="If set, also write total matter P(k) and transfer function files. Default: disabled (only P_cb and redshifts.dat).")
    p.add_argument("--hubble-file", type=str, default="hubble.dat",
                   help="Filename for the Hubble table (two columns: z  H(z) [km/s/Mpc]). Default: hubble.dat")

    args = p.parse_args()

    # P(k) redshift grid (user controlled)
    zs = parse_redshifts(args)
    Nz = len(zs)
    # Background grid (fixed)
    a_bg, z_bg = background_scale_factor_grid()
    os.makedirs(args.outdir, exist_ok=True)

    # Choose k grid & CAMB settings based on units
    h = args.H0 / 100.0
    if args.units == "h":
        k_grid = np.logspace(np.log10(args.kmin), np.log10(args.kmax), args.nk)  # h/Mpc
        kmax_mpc = args.kmax * h  # convert to 1/Mpc for CAMB backend
        hubble_units = True
        k_hunit = True
    else:  # physical units
        k_grid = np.logspace(np.log10(args.kmin), np.log10(args.kmax), args.nk)  # 1/Mpc
        kmax_mpc = args.kmax  # already 1/Mpc
        hubble_units = False
        k_hunit = False

    # CAMB parameters
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=args.H0, ombh2=args.Ob0 * h**2, omch2=args.Ocdm0 * h**2,
                       mnu=args.mnu, nnu=args.Neff, tau=args.tau)
    pars.InitPower.set_params(As=args.As, ns=args.ns)
    pars.NonLinear = model.NonLinear_pk if args.nonlinear else model.NonLinear_none

    # CAMB redshift list: include all non-negative P(k) redshifts plus the maximum background redshift for H(z)
    pk_zs = [float(z) for z in zs if z >= 0.0]
    if len(pk_zs) == 0:
        pk_zs = [0.0]
    z_bg_max = float(np.max(z_bg))  # very high (~9999)
    if z_bg_max not in pk_zs:
        camb_redshifts = sorted(set(pk_zs + [z_bg_max]), reverse=False)
    else:
        camb_redshifts = pk_zs
    pars.set_matter_power(redshifts=camb_redshifts, kmax=kmax_mpc)  # enables transfer + P(k)

    # Run CAMB once
    results = camb.get_results(pars)

    # Interpolator for P_cb,cb with hubble_units=True and k_hunit=True
    PK_cb = results.get_matter_power_interpolator(
        nonlinear=args.nonlinear,
        var1="delta_nonu",  # CDM+baryons
        var2="delta_nonu",
        hubble_units=hubble_units,
        k_hunit=k_hunit
    )
    PK_tot = results.get_matter_power_interpolator(
        nonlinear=args.nonlinear,
        var1="delta_tot",   # CDM+baryons+massive nu
        var2="delta_tot",
        hubble_units=hubble_units,
        k_hunit=k_hunit
    )

    # Matter transfer data (linear). Provides arrays per redshift and k/h grid.
    # CAMB always stores the wavenumbers for transfer data as k/h (h/Mpc).
    mt = results.get_matter_transfer_data()
    for i, z in enumerate(zs):
        # Always write P_cb
        if z >= 0.0:
            Pk_cb_vals = PK_cb.P(z, k_grid)
            pk_cb_fname = os.path.join(args.outdir, f"{args.pk_prefix}_cb_{i:03d}.dat")
            with open(pk_cb_fname, "w") as f:
                for kk, pp in zip(k_grid, Pk_cb_vals):
                    f.write(f"{kk:.8e} {pp:.8e}\n")

        if args.write_optional_files and z >= 0.0:
            # Total matter P(k)
            Pk_tot_vals = PK_tot.P(z, k_grid)
            pk_tot_fname = os.path.join(args.outdir, f"{args.pk_prefix}_tot_{i:03d}.dat")
            with open(pk_tot_fname, "w") as f:
                for kk, pp in zip(k_grid, Pk_tot_vals):
                    f.write(f"{kk:.8e} {pp:.8e}\n")

            # Transfer functions
            zgrid = np.asarray(results.transfer_redshifts)
            iz = int(np.abs(zgrid - z).argmin())
            if abs(zgrid[iz] - z) > 1e-6:
                print(f"[warn] using nearest transfer z={zgrid[iz]:.6f} for requested z={z:.6f}")
            kh_native = mt.transfer_z("k/h", z_index=iz)
            k_native = kh_native * h if args.units == "physical" else kh_native
            T_native = [mt.transfer_z(label, z_index=iz) for label in ["delta_cdm", "delta_baryon", "delta_photon", "delta_neutrino", "delta_nu", "delta_tot"]]
            mask = k_native > 0
            Tk = [np.interp(np.log(k_grid), np.log(k_native[mask]), T[mask]) for T in T_native]
            tf_fname = os.path.join(args.outdir, f"{args.tf_prefix}_{i:03d}.dat")
            with open(tf_fname, "w") as f:
                for kk, tt0, tt1, tt2, tt3, tt4, tt5 in zip(k_grid, Tk[0], Tk[1], Tk[2], Tk[3], Tk[4], Tk[5]):
                    f.write(f"{kk:.8e} {tt0:.8e} {tt1:.8e} {tt2:.8e} {tt3:.8e} {tt4:.8e} {tt5:.8e}\n")
    # P(k) redshift index file (unchanged semantics)
    z_fname = os.path.join(args.outdir, "redshifts.dat")
    with open(z_fname, "w") as f:
        for i, zz in enumerate(zs):
            f.write(f"{i:03d} {zz:.8e}\n")

    # Write Hubble table: two columns z  H(z) [km/s/Mpc]
    hubble_path = os.path.join(args.outdir, args.hubble_file)
    with open(hubble_path, "w") as f:
        for zz in z_bg:
            Hz = results.hubble_parameter(zz)/args.H0  # km/s/Mpc
            f.write(f"{zz:.8e} {Hz:.8e}\n")

    # Optionally write Omega_cb(z) table: two columns z  Omega_cb(z)
    if args.write_optional_files:
        omega_path = os.path.join(args.outdir, "omega_cb.dat")
        with open(omega_path, "w") as f:
            for zz in z_bg:
                Ez = results.hubble_parameter(zz) / args.H0
                omega_cb = (args.Ob0 + args.Ocdm0) * (1.0 + zz) ** 3 / (Ez * Ez)
                f.write(f"{zz:.8e} {omega_cb:.8e}\n")

    extra = []
    if args.write_optional_files:
        extra.append(f"{Nz} T(k) files with prefix '{args.tf_prefix}'")
        extra.append("omega_cb.dat")
    extras = ((" and " + " and ".join(extra)) if extra else "")
    print(
        "Wrote "
        f"{Nz} P(k) files (prefix '{args.pk_prefix}'){extras}, redshifts.dat, "
        f"background H(z) grid ({args.hubble_file}) and {'omega_cb.dat, ' if args.write_optional_files else ''}"
        f"using fixed a-grid (200 points 0.0001->1.51356). Output dir: {os.path.abspath(args.outdir)} (units: {args.units})."
    )

if __name__ == "__main__":
    main()
