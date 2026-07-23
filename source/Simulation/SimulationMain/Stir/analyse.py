#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Christoph Federrath, 2022-2024

import argparse
import numpy as np
import cfpack as cfp
from cfpack import print, stop
import flashlib as fl

# =========================================================
def Gaussian_func(x, mean, sigma):
    ret = 1.0/np.sqrt(2.0*np.pi*sigma**2) * np.exp(-0.5*(x-mean)**2/sigma**2)
    return ret
# =========================================================

# === analyse turbulent field in HDF5 file ===
def analyse(filename):

    print("Analysing file '"+args.inputfile+"'", color="magenta")

    # read data
    gg = fl.FlashGG(filename)
    dirs = ["velx", "vely", "velz"]
    dirs_label = ["$v_x$", "$v_y$", "$v_z$"]
    dat = []
    for dir in dirs:
        dat.append(gg.GetUniformGrid(dir))
    dat = np.array(dat)
    ncmp = len(dirs)

    # statistics
    mean = np.array([dat[d].mean() for d in range(ncmp)])
    std = np.array([dat[d].std() for d in range(ncmp)])
    print("mean velocity = ", mean, highlight=1)
    print("standard deviation of velocity = ", std, highlight=1)
    # PDF plot
    for d in range(ncmp):
        pdfo = cfp.get_pdf(dat[d], bins=min(100,int(0.1*dat[d].size)))
        q = pdfo.bin_center
        pdf = pdfo.pdf
        cfp.plot(x=q, y=pdf, label=dirs_label[d])
        # fit and plot fit
        print("fit for "+dirs[d]+" component:")
        params = {"mean": [-1.0, 0.0, 1.0], "sigma": [1e-6, 1.0, 1e6]}
        fit_result = cfp.fit(Gaussian_func, q, pdf, params=params)
        xfit = np.linspace(q.min(),q.max(),num=500)
        yfit = Gaussian_func(xfit, *fit_result.popt)
        cfp.plot(x=xfit, y=yfit, label="$\mathrm{fit:}\:$"+dirs_label[d])
    # finally, plot everything
    cfp.plot(xlabel='$v$', ylabel='PDF', ylog=True, legend_loc="lower left", save=filename+"_pdf.pdf")

    # Fourier spectra power comparison
    sp = cfp.get_spectrum(dat, ncmp=ncmp)
    # Fourier spectra plot
    k = sp["k"]
    ind = (k > 0) & (sp["P_tot"] > 1e-15)
    cfp.plot(x=k[ind], y=sp["P_tot"][ind], label="total")
    if "P_lgt" in sp: cfp.plot(x=k[ind], y=sp["P_lgt"][ind], label="longitudinal")
    if "P_trv" in sp: cfp.plot(x=k[ind], y=sp["P_trv"][ind], label="transverse")
    cfp.plot(xlabel='$k$', ylabel='$P_v(k)$', xlim=[0.9,1.1*np.max(k)], xlog=True, ylog=True, legend_loc='upper right', save=filename+"_spectrum.pdf")

# =========================================================


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    # read arguments
    parser = argparse.ArgumentParser(description='FLASH post-processing analysis script.')
    parser.add_argument('-i', nargs="?", dest="inputfile", default="Turb_hdf5_plt_cnt_0030", help="input flash file")
    args = parser.parse_args()

    analyse(args.inputfile)
