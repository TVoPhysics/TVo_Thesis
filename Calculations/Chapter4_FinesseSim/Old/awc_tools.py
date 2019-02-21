#!/usr/bin/env python

# This is a set of common tools to use with the general LIGO like files we have.
# Run 'python [filename] -h' to see a list of options.

import sys
from optparse import OptionParser

default_kat_file = "aLIGO_IFO_AWC_tuning.kat"

Lock_Block = "locks"

DOFs  = ["PRCL", "MICH", "CARM", "DARM", "SRCL"]

DOFErrSigs = {"CARM": "REFL_f1_I",
              "DARM": "AS_f2_I",
              "PRCL": "POP_f1_I",
              "MICH": "POP_f2_Q",
              "SRCL": "REFL_f2_I"}

DOFCommands = {}

DOFCommands["CARM"] = """
variable CARM 0
xaxis CARM phi lin -0.001 0.001 10
put* ETMXHR phi $x1
put* ETMXAR phi $x1
put* ETMYHR phi $x1
put* ETMYAR phi $x1"""

DOFCommands["DARM"] = """
variable DARM 0
xaxis DARM phi lin -0.01 0.01 10
put* ETMXHR phi $x1
put* ETMXAR phi $x1
put* ETMYHR phi $mx1
put* ETMYAR phi $mx1"""

DOFCommands["MICH"] = """
variable MICH 0
xaxis MICH phi lin -0.5 0.5 10

put* ITMXHR phi $x1
put* ITMXAR phi $x1
put* ETMXHR phi $x1
put* ETMXAR phi $x1

put* ITMYHR phi $mx1
put* ITMYAR phi $mx1
put* ETMYHR phi $mx1
put* ETMYAR phi $mx1"""

DOFCommands["PRCL"] = """
xaxis* PRMHR phi lin -0.1 0.1 10
put PRMAR phi $x1"""

DOFCommands["SRCL"] = "xaxis* SRMHR phi lin -1 1 10"

Powers = {"DARM": ("Psrc", "Pprc", "asc0", "asf1", "asf2", "src0", "srcf1", "srcf2", "P_DC_OMC"),
          "CARM": ("Px", "Py", "Pprc", "Psrc"),
          "PRCL": ("Pprc", "Psrc", "prc0", "prcf1", "prcf2"),
          "SRCL": ("src0", "srcf1", "srcf2"),
          "MICH": ("asc0", "asf1", "asf2")}


def zero_locks(kat_in, verbose=True, LockIterations = 20):
    kat = kat_in.deepcopy()
    
    # setting noxaxis will override any xaxis commands read in
    kat.noxaxis = True
    kat.verbose = False
    out = kat.run(printerr=1)
    
    # Select lock values from the output
    locks = {x:float(out[x]) for x in out.ylabels if "lock" in x}

    for n in range(0, LockIterations):
        kat.ETMXAR.phi += locks["MICH_lock"] + locks["DARM_lock"] + locks["CARM_lock"]
        kat.ETMXHR.phi += locks["MICH_lock"] + locks["DARM_lock"] + locks["CARM_lock"]
        kat.ETMYAR.phi += locks["CARM_lock"] - locks["MICH_lock"] - locks["DARM_lock"]
        kat.ETMYHR.phi += locks["CARM_lock"] - locks["MICH_lock"] - locks["DARM_lock"]
        kat.ITMXAR.phi += locks["MICH_lock"]
        kat.ITMXHR.phi += locks["MICH_lock"]
        kat.ITMYAR.phi -= locks["MICH_lock"]
        kat.ITMYHR.phi -= locks["MICH_lock"]
        kat.PRMAR.phi  += locks["PRCL_lock"]
        kat.PRMHR.phi  += locks["PRCL_lock"]
        kat.SRMHR.phi  += locks["SRCL_lock"]
        
        # Get the new lock values, should be closer to 0...
        out = kat.run(printerr=1, cmd_args=["-cr=on"])
        locks = {x:float(out[x]) for x in out.ylabels if "lock" in x}

    if verbose:
        from pprint import pprint # used to print the lock dictionary, not really needed
        
        print("")
        print ("Final lock values:")
        
        # print ditionary of values one per line
        pprint(locks, width=1) 
        
        print ("Final error signal values:")
        
        for DOF in DOFs:
            print("%s : %s = %.15g" % (DOF, DOFErrSigs[DOF], out[DOFErrSigs[DOF]]))

        print ("Final tunings:")
        print ("")
        print ("const phi_SRM %.15g" % kat.SRMHR.phi)
        print ("const phi_PRM %.15g" % kat.PRMHR.phi)
        print ("const phi_ITMX %.15g" % kat.ITMXHR.phi)
        print ("const phi_ITMY %.15g" % kat.ITMYHR.phi)
        print ("const phi_ETMX %.15g" % kat.ETMXHR.phi)
        print ("const phi_ETMY %.15g" % kat.ETMYHR.phi)
        print ("")


def plot_DOFs(_kat, savefig=None, steps=10):
    import pykat
    import pylab
    import scipy
    
    import matplotlib.gridspec as gridspec
    from scipy.interpolate import interp1d
    from matplotlib.ticker import EngFormatter
    
    _kat.verbose = False
    
    _kat.removeBlock(Lock_Block, False)

    pylab.figure(figsize=(10,8))
    gs = gridspec.GridSpec(3, 2)
    ax = [None, None, None, None, None]

    ax[0] = pylab.subplot(gs[0, 0])
    ax[1] = pylab.subplot(gs[0, 1])
    ax[2] = pylab.subplot(gs[1, 0])
    ax[3] = pylab.subplot(gs[1, 1])
    ax[4] = pylab.subplot(gs[2, :])

    subplot = 0

    for DOF in DOFs:
        kat = _kat.deepcopy()
        kat.parseCommands(DOFCommands[DOF])
        kat.xaxis.steps = steps
    
        out = kat.run(printerr=1)
    
        f1 = interp1d(out[DOFErrSigs[DOF]], out.x)
        f2 = interp1d(out.x, out[DOFErrSigs[DOF]])
    
        try:
            title = "0-Crossing = %.3g\nDeg, offset at 0 = %.3g" % (f1(0), f2(0))
        except:
            title = "Interpolation error"
    
        ax[subplot].ticklabel_format(style="sci", scilimits=(1,2))
        ax[subplot].set_title(title)
        ax[subplot].plot(out.x, out[DOFErrSigs[DOF]])
        ax[subplot].set_xlim(out.x.min(), out.x.max())
        ax[subplot].set_xlabel(DOF + " [Deg]")
        ax[subplot].set_ylabel(DOFErrSigs[DOF])
        ax[subplot].grid()
    
        subplot += 1

    _kat.noxaxis = True
    out = _kat.run()

    detectors = [kat.Px, kat.Py, kat.Psrc, kat.Pprc, kat.P_DC_OMC]

    print("\n--- Detectors")

    for det in detectors:
        print("%s = %g W" %(det, out[det]))
    
    pylab.tight_layout()
    
    if savefig is not None:
        fname = str(savefig)
        print("Saving plots to " + fname)
        pylab.savefig(fname)
        
    pylab.show()



def length_sensing_matrix(basekat):

    import pykat
    import numpy as np
    import pandas
    import pylab

    from pandas import DataFrame
    from pykat import beam_param
    from scipy.interpolate import interp1d

    basekat.verbose = False
    basekat.removeBlock(Lock_Block, False)
    basekat.yaxis = "re:im"

    ErrSigNames = [DOFErrSigs[dof] for dof in DOFs]

    SensM = DataFrame(columns=ErrSigNames, index=DOFs)
    
    # Now loop over each DOF 
    for DOF in DOFs:
        print("Computing %s..." % DOF)
    
        kat = basekat.deepcopy()      # Create a copy we can manipulate
        kat.parseCommands(DOFCommands[DOF]) # Load in the DOF specific code
    
        out = kat.run(printerr=1)
    
        legends = []
    
        # Compute the (cross-)coupling between other DOFs
        for errsig in ErrSigNames:
            # We want to compute the gradient of the error signal
            x = out.x 
            y = out[errsig].real
        
            dx   = x[1] - x[0]
            dydx = interp1d(x, np.gradient(y, dx))
        
            # Seems the pandas floating point precision setting is just a hint
            # and it somtimes ignores it. So here we force 4 d.p. for printing
            # the matrix
            SensM[errsig][DOF] = dydx(0)
    
    print(SensM)

    G = DataFrame(-1.0/SensM.as_matrix(), columns=ErrSigNames, index=DOFs)
    
    print("\nNew lock commands to include in file (Note: Lock accuracy set to a default value!)\n")
    
    for s in zip(ErrSigNames, DOFs):
        print("lock {DOF}_lock ${DOF}_err {Gain} 10u".format(DOF=s[1], Gain=G[s[0]][s[1]]))
    
    print("")


def plot_powers(_kat, savefig=None, steps=100):
    import pykat
    import pylab
    import scipy
    import matplotlib.gridspec as gridspec
    from scipy.interpolate import interp1d
    from matplotlib.ticker import EngFormatter
    
    pylab.figure(figsize=(15,10))

    _kat.verbose = False
    _kat.removeBlock(Lock_Block, False)
    _kat.maxtem = 0

    gs = gridspec.GridSpec(3, 2)
    ax = [None, None, None, None, None]

    ax[0] = pylab.subplot(gs[0, 0])
    ax[1] = pylab.subplot(gs[0, 1])
    ax[2] = pylab.subplot(gs[1, 0])
    ax[3] = pylab.subplot(gs[1, 1])
    ax[4] = pylab.subplot(gs[2, :])

    subplot = 0

    for DOF in DOFs:
        kat = _kat.deepcopy()
        kat.parseCommands(DOFCommands[DOF])
        kat.xaxis.steps = steps
        
        kat.xaxis.limits[0] *= 1
        kat.xaxis.limits[1] *= 1
        
        out = kat.run(printerr=1)
    
        ax[subplot].ticklabel_format(style="sci", scilimits=(1,2))
        ax[subplot].set_title(DOF)
    
        for det in Powers[DOF]:
            ax[subplot].semilogy(out.x, out[det], label=det)

        ax[subplot].legend(ncol=2, fontsize="small")
        ax[subplot].set_xlim(out.x.min(), out.x.max())
        ax[subplot].set_xlabel(DOF + " [Deg]")
        ax[subplot].grid()
    
        subplot += 1

    _kat.noxaxis = True
    out = _kat.run()
    
    pylab.tight_layout()
    
    if savefig is not None:
        fname = str(savefig)
        print("Saving plots to " + fname)
        pylab.savefig(fname)
    
    pylab.show()
    
    
    
def dc_offset(kat, dc_off, verbose=False):
    
    import pykat
    import pylab

    from math import pi
    from pprint import pprint # used to print the lock dictionary, not really needed
    
    kat.verbose = False

    # Firstly disable all the detectors
    for x in kat.detectors: kat.detectors[x].enabled = False

    # enable those that the locks need
    for pd in [kat.detectors[pd] for pd in DOFErrSigs.values()]:
        pd.enabled = True

    # First we zero all the locks
    if verbose: print("Zeroing locks...")
    zero_locks(kat)

    # extra factor of 0.5 here as DARM is (ETMX-ETMY)/2?
    phi_off = dc_off * (2*pi / 1064e-9) * (180 / pi) / 2.0

    # make a copy of the zero'd kat file to compute DC offset levels
    _kat = kat.deepcopy()

    _kat.DARM_lock.remove()
    _kat.ETMX_lock.remove()
    _kat.ETMY_lock.remove()

    # as we are already putting to ETMX and ETMY we need to 
    # xaxis a dummy variable and put that into the 
    _kat.parseCommands("""
    func ETMY_lock = $CARM_lock - $MICH_lock - $x1
    func ETMX_lock = $CARM_lock + $MICH_lock + $x1

    variable darm 0
    xaxis darm re lin 0 {phi_off} 100
    """.format(phi_off=phi_off) )

    _kat.noxaxis = False
    _kat.verbose = False

    out = _kat.run(printerr=1)
    AS_f2_I_offset = float(out["AS_f2_I"][-1])

    if verbose: 
        print("")
        print("AS_f2_I at %g offset = %g" % (dc_off, AS_f2_I_offset))
        print("")

    # Now go back to the original kat file and
    # remove the original error signal as now we have to apply our offset
    # in a func
    kat.remove("set DARM_err AS_f2_I re")

    # As the ordering of locks and functions is crucial
    # we need to remove all the previous locks depending on DARM
    # and readd them in specific order
    DARM_lock = kat.DARM_lock.getFinesseText()
    ETMX_lock = kat.ETMX_lock.getFinesseText()
    ETMY_lock = kat.ETMY_lock.getFinesseText()

    kat.DARM_lock.remove()
    kat.ETMX_lock.remove()
    kat.ETMY_lock.remove()

    # Now add a new varriable to read the AS_f2_P signal and
    # a function to apply our offset.
    # As we are feeding a function into another function we have to 
    # make sure that the insert function is before the other in the code
    DARM_offset_cmds = """
    set AS_f2_I_re AS_f2_I re
    func DARM_err = $AS_f2_I_re - {offset}
    {DARM_lock}
    {ETMX_lock}
    {ETMY_lock}
    """.format(offset=AS_f2_I_offset, DARM_lock=DARM_lock, ETMX_lock=ETMX_lock, ETMY_lock=ETMY_lock)

    kat.parseCommands(DARM_offset_cmds, addToBlock="locks")

    if verbose: print("Zeroing locks again...")
    zero_locks(kat, verbose=True)

    kat.P_DC_OMC.enabled = True
    kat.noxaxis = True
    kat.removeBlock(Lock_Block)
    out = kat.run()

    if verbose: 
        print("")
        print("DC power with offset = %s W" % float(out[kat.P_DC_OMC]))
        print("")
        print("DARM error signal with offset commands to include in kat:")
        print(DARM_offset_cmds)    



if __name__ == "__main__":
    
    parser = OptionParser()
    
    parser.add_option("-k", "--kat",
                      dest="kat_file", help="Sets kat file to run. Current default is '%s'" % default_kat_file, default=default_kat_file)

    parser.add_option("-s", "--sensing-matrix", dest="length_sensing_matrix", action="store_true",
                      help="Plot powers in the various cavities", default=False)
                      
    parser.add_option("-z", "--zero", action="store_true",
                      dest="zero_file", help="Zero locks of kat file", default=False)

    parser.add_option("-d", "--plot-dof", dest="plot_dof", action="store_true",
                      help="Plots error signals over DOFs", default=False)

    parser.add_option("-p", "--plot-powers", dest="plot_powers", action="store_true",
                      help="Plot powers in the various cavities", default=False)
                      
    parser.add_option("-o", "--dc-offset", dest="dc_offset",
                      help="Computes values for a DC offset", default=None)
                      
    (options, args) = parser.parse_args()

    print("Using kat file " + options.kat_file)
    
    if options.zero_file:
        import pykat
        kat = pykat.finesse.kat()
        kat.verbose = False
        kat.loadKatFile(options.kat_file)
        zero_locks(kat)
        
    if options.plot_dof:
        import pykat
        kat = pykat.finesse.kat()
        kat.verbose = False
        kat.loadKatFile(options.kat_file)
        plot_DOFs(kat)
    
    if options.plot_powers:
        import pykat
        kat = pykat.finesse.kat()
        kat.verbose = False
        kat.loadKatFile(options.kat_file)
        plot_powers(kat)
    
    if options.length_sensing_matrix:
        import pykat
        kat = pykat.finesse.kat()
        kat.verbose = False
        kat.loadKatFile(options.kat_file)
        length_sensing_matrix(kat)
    
    if options.dc_offset is not None:
        import pykat
        kat = pykat.finesse.kat()
        kat.verbose = False
        kat.loadKatFile(options.kat_file)
        dc_offset(kat, float(options.dc_offset), verbose=True)
        
    
    
    