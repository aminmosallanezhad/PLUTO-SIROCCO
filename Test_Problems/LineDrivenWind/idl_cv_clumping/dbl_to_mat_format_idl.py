import os
import sys
from numpy import *
from matplotlib.pyplot import *
import scipy.io as io
import pyPLUTO.pload as pp

entryRange = 1000
i = 0


while i <= entryRange:
    D = pp.pload(i, datatype='dbl')
    io.savemat('data.{}.mat'.format(i),{'M_max1': D.M_max1,'M_max2': D.M_max2, 'T': D.T, 'T_r': D.T_r, 'XI': D.XI, 'bc': D.bc, 'bc_pre': D.bc_pre, 'cc': D.cc, 'cc_pre': D.cc_pre, 'ch': D.ch, 'ch_pre': D.ch_pre, 'dv_ds' :D.dv_ds, 'gr': D.gr, 'gt': D.gt, 'lc': D.lc, 'lc_pre': D.lc_pre, 'ne': D.ne, 'nh': D.nh, 'rho' :D.rho, 'prs' :D.prs, 't' :D.t, 'tr1': D.tr1, 'vx1': D.vx1, 'vx2': D.vx2, 'vx3': D.vx3, 'xh' :D.xh, 'xh_pre' :D.xh_pre, 'x1' :D.x1, 'x2':D.x2 , 'dx1' :D.dx1, 'dx2':D.dx2 , 'x1r' :D.x1r, 'x2r':D.x2r})
    i = i+1
