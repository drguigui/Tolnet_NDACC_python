#!/usr/bin/env python3

from pylab import *




def NDAAC_ResolDf(vZsampling, vCoef, vHfinp):
    """ Computes the resolution using the FWHM of the impulse response

    ... (add Thierry data)


    Returns Status, 
    Returns dzfwhm: the resolution in the same unit of vZsampling
    Returns irout: the impulse response
    """
    ncoef = len(vCoef)
    nf = len(vHfinp)
    status = 0
    missval = -99

    hfinp = vHfinp.copy()

    hf = zeros((nf)) 
    hfout = hf.copy()
    dzcutf = missval
    fc = missval
    fcinp = missval

    if not (ncoef % 2): # the number of coefs is even, we reject
        print("Error, The number of coefficients is even")
        return status, 0, irout

    if ncoef == 1 and vCoef[0] != 1: # We reject when one coeff different from 1
        return status, 0, irout

    if ncoef == 1 and vCoef[0] == 1: # We reject when one coeff different from 1
        status = 1
    else: # at least 3 coefs (number odd) => we calculate the transfer functions 
        nc = (ncoef - 1) // 2
        # we check if the coefficients are symmetric (status 1) or antisymmetric (status 2)
        status = -1
        for ic in range(1, nc + 1): # Python is like IDL, subscript starts at 0
            if abs(vCoef[nc - ic]) < 1E-10 and abs(vCoef[nc + ic]) < 1E-10:
                print("0th coeff")
            elif abs(vCoef[nc - ic] - vCoef[nc + ic]) < 1E-3 * abs(vCoef[nc + ic]):
                if status == 2:
                    status = 0
                    print("error df, we have a mix of symmetric and antisymmetric")
                    return status, 0, irout
                status = 1
            elif (abs(vCoef[nc - ic] + vCoef[nc + ic]) < 1E-3 * abs(vCoef[nc + ic]) and vCoef[nc] == 0):
                if status == 1:
                    status = 0
                    print("error df, we have a mix of symmetric and antisymmetric")
                    return status, 0, irout
                status = 2
            else:
                    status = 0
                    print("error df , we have a mix of symmetric and antisymmetric")
                    return status, 0, irout


    f = 0.5 / nf + 0.5 / nf * array(list(range(nf)))

    for jf in range(nf - 1):
        if (hfinp[jf + 1] - 0.5) * (hfinp[jf] - 0.5) <= 0:
            if abs(hfinp[jf] - 0.5) < 1E-9: #TODO !!!! -9
                fcinp = f[jf]
            elif  abs(hfinp[jf + 1] - 0.5) < 1E-9: #TODO !!!! -9: 
                fcinp = f[jf + 1]
            else:
                weight = abs(hfinp[jf] - 0.5) / abs(hfinp[jf + 1] - hfinp[jf])
                fcinp = f[jf] + weight * abs(f[jf+1] - f[jf])
            break

    if sum(hfinp) == nf and hfinp.std() == 0:
        fcinp = 0.5


    hfcomp = zeros((nf), dtype=complex)
    for jf in range(nf):
        summ = 0 + 0.j
        for ic in range(-nc, nc + 1):
            summ += vCoef[nc - ic] * (cos(2 * pi * f[jf] * ic) - 1.j * sin(2 * pi * f[jf] * ic))
        hfcomp[jf] = summ
    hf = hfcomp.copy()
    hfout = hf * hfinp
    if status == 2:
        #print hfout #, complex(hfout)

        hfout = imag(hfout) / (2 * pi * f)
    else:
        hfout = real(hfout)
    
    hfc = 0.5
    db = -6
    #hfc = 10 ** (db / 20.)

    for jf in range(nf - 1):
        if (hfout[jf + 1] - hfc) * (hfout[jf] - hfc) <= 0:
            if abs(hfout[jf] - hfc) < 1E-9: 
                fc = f[jf]
            elif  abs(hfout[jf + 1] - hfc) < 1E-9:
                fc = f[jf + 1]
            else:
                weight = abs(hfout[jf] - hfc) / abs(hfout[jf + 1] - hfout[jf])
                fc = f[jf] + weight * abs(f[jf+1] - f[jf])
            break

    if ncoef == 1:
        fc = fcinp
    if fc > 0:
        dzcutf = vZsampling / (2 * fc)
    return status, dzcutf, hfout


if '__main__' == __name__:
    impulse = ones((2400))
   # status, dz, impulse =  NDAAC_ResolDf(1, array([1/5,1/5,1/5,1 /5.,1/5]), impulse)
    
    #impulse = ones((2400))
   # print NDAAC_ResolDf(1, 1/9*array([1,1,1,1,1,1,1,1,1]), impulse)

  #  print  NDAAC_ResolDf(1, array([-3,-2,-1,0,1,2,3]) / 28., impulse)
    #print  NDAAC_ResolDf(1, ones((25)) / 25, impulse)
    #status, dz, impulse =  NDAAC_ResolDf(7.5, ones((25)) / 25, impulse)
    #print NDAAC_ResolDf(1, 1/2. * array([-1,0,1]), impulse)
    impulse = ones((2400))
    result =  NDAAC_ResolDf(7.5 * 25, 1/2. * array([-1,0,1]), impulse)
    print(result[1])
   # print 7.5 * 25



