#!/usr/bin/env python3
"""
Compute both P_cb(k,z) and T_cb(k,z) with CAMB and write one pair of files per redshift.

Outputs per redshift index i:
  pk_cb_{i:03d}.txt  -> columns: k_h [h/Mpc], P_cb [(Mpc/h)^3]
  tf_cb_{i:03d}.txt  -> columns: k_h [h/Mpc], T_cb [dimensionless]

Notes:
- Power spectrum uses CAMB's interpolator with var1=var2='delta_nonu' (CDM+baryons).
- Transfer function is taken from results.get_matter_transfer_data() for 'delta_nonu'.
- Transfer is interpolated onto the same k_h grid used for P(k) output.
"""

import argparse
import os
import numpy as np

import camb
from camb import model

def parse_redshifts(args) -> np.ndarray:
    if args.zs is not None:
        zs = [float(x) for x in args.zs.split(",") if x.strip() != ""]
        return np.array(zs, dtype=float)
    elif args.nz is not None and args.zmin is not None and args.zmax is not None:
        return np.linspace(args.zmin, args.zmax, args.nz, dtype=float)[::-1]
    else:
        return np.array([0.0, 0.5, 1.0], dtype=float)[::-1]

def main():
    p = argparse.ArgumentParser(description="Compute CDM+baryons P(k) and transfer T(k) with CAMB.")
    # Redshift controls
    p.add_argument("--zs", type=str, default=None,
                   help="Comma-separated list of redshifts, e.g., '0,0.5,1'.")
    p.add_argument("--zmin", type=float, default=None, help="Minimum redshift (used with --zmax/--nz).")
    p.add_argument("--zmax", type=float, default=None, help="Maximum redshift (used with --zmin/--nz).")
    p.add_argument("--nz", type=int, default=None, help="Number of redshifts (used with --zmin/--zmax).")

    # k-grid in h/Mpc
    p.add_argument("--kmin-h", type=float, default=1e-4, help="Minimum k in h/Mpc (default: 1e-4).")
    p.add_argument("--kmax-h", type=float, default=10.0, help="Maximum k in h/Mpc (default: 10).")
    p.add_argument("--nk", type=int, default=512, help="Number of k points (log-spaced; default: 512).")

    # Cosmology (Planck-like defaults)
    p.add_argument("--H0", type=float, default=67.66, help="H0 in km/s/Mpc (default: 67.66).")
    p.add_argument("--ombh2", type=float, default=0.02237, help="Omega_b h^2 (default: 0.02237).")
    p.add_argument("--omch2", type=float, default=0.1200, help="Omega_c h^2 (default: 0.1200).")
    p.add_argument("--ns", type=float, default=0.9649, help="Scalar spectral index n_s (default: 0.9649).")
    p.add_argument("--As", type=float, default=2.1e-9, help="Primordial amplitude A_s (default: 2.1e-9).")
    p.add_argument("--mnu", type=float, default=0.06, help="Sum of neutrino masses in eV (default: 0.06).")
    p.add_argument("--Neff", type=float, default=3.044, help="Effective N_eff for massless species (default: 3.044).")
    p.add_argument("--tau", type=float, default=0.054, help="Optical depth tau (default: 0.054).")

    # Linear / non-linear choice for P(k)
    p.add_argument("--nonlinear", action="store_true",
                   help="If set, include non-linear corrections for matter P(k). Linear by default.")

    # Output controls
    p.add_argument("--outdir", type=str, default=".", help="Output directory (default: current).")
    p.add_argument("--pk-prefix", type=str, default="pk_cb_", help="P(k) filename prefix (default: 'pk_cb_').")
    p.add_argument("--tf-prefix", type=str, default="tf_cb_", help="T(k) filename prefix (default: 'tf_cb_').")

    args = p.parse_args()

    zs = parse_redshifts(args)
    Nz = len(zs)
    os.makedirs(args.outdir, exist_ok=True)

    # k grid in h/Mpc for evaluation and output
    kh = np.logspace(np.log10(args.kmin_h), np.log10(args.kmax_h), args.nk)

    # CAMB parameters
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=args.H0, ombh2=args.ombh2, omch2=args.omch2,
                       mnu=args.mnu, nnu=args.Neff, tau=args.tau)
    pars.InitPower.set_params(As=args.As, ns=args.ns)
    pars.NonLinear = model.NonLinear_pk if args.nonlinear else model.NonLinear_none

    # CAMB wants kmax in 1/Mpc (not h/Mpc); convert using h = H0/100
    h = args.H0 / 100.0
    kmax_mpc = args.kmax_h * h

    # Ensure CAMB computes transfers and matter power at requested redshifts and k-range
    pars.set_matter_power(redshifts=list(zs), kmax=kmax_mpc)  # computes transfer data and enables P(k)

    # Run CAMB once
    results = camb.get_results(pars)

    # Interpolator for P_cb,cb with hubble_units=True and k_hunit=True
    PK = results.get_matter_power_interpolator(
        nonlinear=args.nonlinear,
        var1="delta_nonu",  # CDM+baryons
        var2="delta_nonu",
        hubble_units=True,  # P in (Mpc/h)^3
        k_hunit=True        # k in h/Mpc
    )

    # Matter transfer data (linear). Provides arrays per redshift and k/h grid.
    mt = results.get_matter_transfer_data()
    # Mapping from requested P(k) redshifts to internal transfer redshift indices
    pk_to_transfer_idx = results.PK_redshifts_index

    for i, z in enumerate(zs):
        # ----- Power spectrum P_cb(k,z) on chosen k_h grid -----
        Pk = PK.P(z, kh)
        pk_fname = os.path.join(args.outdir, f"{args.pk_prefix}{i:03d}.dat")
        with open(pk_fname, "w") as f:
            #f.write(f"# P_cb,cb(k,z) from CAMB\n")
            #f.write(f"# cosmology: H0={args.H0}, ombh2={args.ombh2}, omch2={args.omch2}, mnu={args.mnu}, Neff={args.Neff}, ns={args.ns}, As={args.As}\n")
            #f.write(f"# nonlinear={'True' if args.nonlinear else 'False'}\n")
            #f.write(f"# z = {z}\n")
            #f.write(f"# k_h [h/Mpc]    P_cbcb [(Mpc/h)^3]\n")
            for kk, pp in zip(kh, Pk):
                f.write(f"{kk:.8e} {pp:.8e}\n")

        # ----- Transfer function T_cb(k,z) on the same k_h grid -----
        # Map requested z to nearest internal transfer redshift
        zgrid = np.asarray(results.transfer_redshifts)          # CAMB's transfer z's
        iz = int(np.abs(zgrid - z).argmin())

        # Optional sanity check
        if abs(zgrid[iz] - z) > 1e-6:
            print(f"[warn] using nearest transfer z={zgrid[iz]:.6f} for requested z={z:.6f}")

        kh_native = mt.transfer_z("k/h", z_index=iz)
        Tcb_native = mt.transfer_z("delta_nonu", z_index=iz)     # CDM+baryons transfer

        # Interpolate onto your chosen kh grid (log-k interpolation)
        mask = kh_native > 0
        Tcb = np.interp(np.log(kh), np.log(kh_native[mask]), Tcb_native[mask])

        tf_fname = os.path.join(args.outdir, f"{args.tf_prefix}{i:03d}.dat")
        with open(tf_fname, "w") as f:
            #f.write(f"# T_cb(k,z) from CAMB (CDM+baryons)\n")
            #f.write(f"# cosmology: H0={args.H0}, ombh2={args.ombh2}, omch2={args.omch2}, mnu={args.mnu}, Neff={args.Neff}, ns={args.ns}, As={args.As}\n")
            #f.write(f"# linear transfer; requested z = {z} ; internal z = {zgrid[iz]}\n")
            #f.write(f"# k_h [h/Mpc]    T_cb [dimensionless]\n")
            for kk, tt in zip(kh, Tcb):
                f.write(f"{kk:.8e} {tt:.8e} {tt:.8e} {tt:.8e} {tt:.8e} {tt:.8e} {tt:.8e}\n")
    z_fname = os.path.join(args.outdir, f"redshifts.dat")
    with open(z_fname, "w") as f:
        for i, zz in enumerate(zs):
            f.write(f"{i:03d} {zz:.8e}\n")

    print(f"Wrote {Nz} P(k) files with prefix '{args.pk_prefix}' and {Nz} T(k) files with prefix '{args.tf_prefix}' to {os.path.abspath(args.outdir)}.")

if __name__ == "__main__":
    main()
