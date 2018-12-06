#!/usr/bin/env python2

""" posforce2beadheat.py

Reads positions and momenta from an i-PI run and computes the bead heat component.

Assumes the input files are in xyz format and atomic units, with
prefix.pos_*.xyz and prefix.mom_*.xyz naming scheme.

Syntax:
   posforce2beadheat.py prefix temperature[K]
"""


import numpy as np
import sys
import glob
from ipi.utils.io import read_file
from ipi.engine.beads import Beads
from ipi.utils.depend import dstrip
from ipi.utils.units import unit_to_internal, Constants
from ipi.utils.messages import verbosity
verbosity.level = "low"


def main(prefix, temp):

    T = float(temp)
    fns_pos = sorted(glob.glob(prefix + ".pos*"))
    fns_mom = sorted(glob.glob(prefix + ".mom*"))
    fn_out_heat = prefix + ".heat.xyz"

    # check that we found the same number of positions and forces files
    nbeads = len(fns_pos)
    if nbeads != len(fns_mom):
        print fns_pos
        print fns_mom
        raise ValueError("Mismatch between number of input files for momenta and positions.")

    # print some information
    print 'temperature = {:f} K'.format(T)
    print
    print 'number of beads = {:d}'.format(nbeads)
    print
    print 'positions and forces file names:'
    for fn_pos, fn_mom in zip(fns_pos, fns_mom):
        print '{:s}   {:s}'.format(fn_pos, fn_mom)
    print
    print 'output file names:'
    print fn_out_heat
    print

    temp = unit_to_internal("energy", "kelvin", T)

    # Calculate w_P, which just depends on the temperature and nbeads

    wp=temp/nbeads

    # open input and output files
    ipos = [open(fn, "r") for fn in fns_pos]
    imom = [open(fn, "r") for fn in fns_mom]
    iheat = open(fn_out_heat, "w")

    natoms = 0
    ifr = 0
    while True:

        # print progress
        if ifr % 100 == 0:
            print '\rProcessing frame {:d}'.format(ifr),
            sys.stdout.flush()

        # load one frame
        try:
            for i in range(nbeads):
                ret = read_file("xyz", ipos[i])
                pos = ret["atoms"]
                ret = read_file("xyz", imom[i])
                mom = ret["atoms"]
                if natoms == 0:
                    natoms = pos.natoms
                    beads = Beads(natoms, nbeads)
                    moms = Beads(natoms, nbeads)
                    heat = np.zeros((natoms, 3), float)
                beads[i].q = pos.q
                moms[i].q = mom.q
        except EOFError:
            # finished reading files
            break

        # calculate kinetic energies
        q = dstrip(beads.q)
        p = dstrip(moms.q)
        heat[:] = 0.0
        for j in range(nbeads):
            for i in range(natoms):
                if j==0:
                   jminus=nbeads-1
                   jplus=j+1
                if j==(nbeads-1):
                   jminus=j-1
                   jplus=0
                else:
                   jplus=j+1
                   jminus=j-1
                heat[i, 0] += ((q[jplus, i * 3 + 0] - q[j, i * 3 + 0]) * (-q[jplus, i * 3 + 0] + q[j, i * 3 + 0]) * p[j, i * 3 + 0] + (q[jminus, i * 3 + 0] - q[j, i * 3 + 0]) * (-q[jminus, i * 3 + 0] + q[j, i * 3 + 0]) * p[j, i * 3 + 0]  )
                heat[i, 1] += ((q[jplus, i * 3 + 1] - q[j, i * 3 + 1]) * (-q[jplus, i * 3 + 1] + q[j, i * 3 + 1]) * p[j, i * 3 + 1] + (q[jminus, i * 3 + 1] - q[j, i * 3 + 1]) * (-q[jminus, i * 3 + 1] + q[j, i * 3 + 1]) * p[j, i * 3 + 1]  )
                heat[i, 2] += ((q[jplus, i * 3 + 0] - q[j, i * 3 + 2]) * (-q[jplus, i * 3 + 2] + q[j, i * 3 + 2]) * p[j, i * 3 + 2] + (q[jminus, i * 3 + 2] - q[j, i * 3 + 2]) * (-q[jminus, i * 3 + 2] + q[j, i * 3 + 2]) * p[j, i * 3 + 2]  )
        heat *= wp**2 / nbeads
#        kcv[:, 0:3] += 0.5 * Constants.kb * temp
#        kcv[:, 3:6] *= 0.5
        fullheat = np.sum(heat, axis=0)
        # write output
        iheat.write("%d\n# Full bead heat flux [a.u.] (x, y, z) = %12.5e %12.5e %12.5e\n" % (natoms, fullheat[0], fullheat[1], fullheat[2]))
        for i in range(natoms):
            iheat.write("%8s %12.5e %12.5e %12.5e\n" % (pos.names[i], heat[i, 0], heat[i, 1], heat[i, 2]))

        ifr += 1

    print '\rProcessed {:d} frames.'.format(ifr)

    iheat.close()

if __name__ == '__main__':
    main(*sys.argv[1:])
