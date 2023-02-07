import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys


fnames = sys.argv[1:]

print("Opening {} HDF5 files...".format(len(fnames)))

chis = []
gtildes = []
# gs = []

xs = []
ts = []

for fname in fnames:
    with h5py.File(fname, "r") as hf:
        ts.append(hf["real scalars"][0][1])
        chis.append(hf["z400"][()])
        gtildes.append(hf["z401"][()])
        xs.append(hf["bounding box"][()])

ts = np.array(ts)
print(ts)

xlo = xs[0][0, 0, 0]
xhi = xs[0][0, 0, 1]

xf = np.linspace(xlo, xhi, chis[0].shape[3]+1)
x = 0.5*(xf[:-1]+xf[1:])
dx = np.diff(x)[0]

nx = x.size

plt.figure()
for i in range(len(chis)):
    chi = chis[i][0, 0, 0, :]
    gtilde = gtildes[i][0, 0, 0, :]

    plt.plot(x, gtilde/chi)


plt.ylabel(r"$\chi$")
plt.xlabel(r"$x$")
plt.savefig("gauge_wave.pdf")

plt.figure()
res = (gtildes[1][0, 0, 0, :]/chis[1][0, 0, 0, :] -
       gtildes[0][0, 0, 0, :]/chis[0][0, 0, 0, :])

plt.plot(x, res)

plt.savefig("gauge_wave_residual.pdf")

nx = chis[0].shape[3]
amplitude = np.max(gtildes[0][0, 0, 0, :]/chis[0][0, 0, 0, :])
L2norm = np.sqrt(np.sum((res/amplitude)**2)) / nx
print("residual norm =", L2norm)
try:
    assert(L2norm < 1e-7)
except:
    print("ERROR: the residual is too large")
else:
    print("Success!")
