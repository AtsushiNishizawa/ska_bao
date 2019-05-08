import pyfits as pf
import numpy as np
import os, sys


def writefits(x_gas, x_dmp, v_gas, v_dmp, nha, z, fits):
    Np = x_gas.shape[1]
    fmt3 = "{0:d}E".format(Np)
    fmt1 = "E"
    extheader = pf.Header()
    extheader["redshift"] = "{0:.3F}".format(z)
    
    coldefs = pf.ColDefs([
            pf.Column( name="x_gas", format=fmt3, array=x_gas ),
            pf.Column( name="x_dmp", format=fmt3, array=x_dmp ),
            pf.Column( name="v_gas", format=fmt3, array=v_gas ),
            pf.Column( name="v_dmp", format=fmt3, array=v_dmp ),
            pf.Column( name="nha"  , format=fmt1, array=nha   )
            ])
    tbhdu = pf.BinTableHDU.from_columns(coldefs, header=extheader)
    tbhdu.writeto(fits, clobber=True)
    print "write to ", fits, Np, "objects"
