#!/usr/bin/env python3

from pylab import *




def NDAAC_ResolIr(vZsampling, vCoef, vIRinp):
    """ Computes the resolution using the FWHM of the impulse response

    ... (add Thierry data)


    Returns Status, 
    Returns dzfwhm: the resolution in the same unit of vZsampling
    Returns irout: the impulse response
    """
    ncoef = len(vCoef)
    nk = len(vIRinp)
    status = 0
    missval = -99
    irout = vIRinp.copy() # we copy the values of vIRinp in irout
    irconvol = vIRinp.copy() # we copy the values of vIRinp in irout
    ir = zeros(nk)  

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
                    print("error ir, we have a mix of symmetric and antisymmetric")
                    return status, 0, irout
                status = 1
            elif (abs(vCoef[nc - ic] + vCoef[nc + ic]) < 1E-3 * abs(vCoef[nc + ic]) and vCoef[nc] == 0):
                if status == 1:
                    status = 0
                    print("error ir, we have a mix of symmetric and antisymmetric")
                    return status, 0, irout
                status = 2
            else:
                    status = 0
                    print("error ir, we have a mix of symmetric and antisymmetric")
                    return status, 0, irout


#    print "0", irconvol
    if sum(irconvol) == 0 and irconvol.std() == 0:
        irconvol[nk // 2] = 1

    if status == 2:
        irtemp = zeros((len(irconvol)))
        irtemp[0] =irconvol[0]
        for k in range(1, nk):
            irtemp[k] = sum(irconvol[0:k + 1])   
        irtemp = irtemp / max(irtemp)
        irconvol = irtemp

#    print irconvol[1190:1210]

    #print "1", ir
    for j in range(nc, len(irconvol) - nc ):
      # ir[j + nc] = sum(array(irconvol[j-nc :j+nc + 1]) * array(vCoef))
       ir[j] = sum(array(irconvol[j-nc :j+nc + 1]) * array(vCoef)) # We have to do the convolution like this
    
    #print "2", ir





    if max(ir) < 0:
        irout = abs(ir)
    else:
        irout = ir
    maxirout = max(irout)
    irout /= maxirout

#    print irout
    db = -6
    #irc = 10 ** (db / 20.)
    irc = 0.5
#    print irc
    
    k1 = 0
    x1 = 0
    x2 = 0
    
    for k in range(nk - 1):
        if (irout[k + 1] - irc) * (irout[k] - irc) <= 0:
            if abs(irout[k] - irc) < 1E-9:
                k1 = k
                x1 = k1
            elif  abs(irout[k + 1] - irc) < 1E-9:
                k1 = k + 1
                x1 = k1
            else:
                weight = abs(irout[k] - irc) / abs(irout[k + 1] - irout[k])
                k1 = k
                x1 = k + weight
            break

    # print "k1, x1", k1, x1
    for k in range(k1 + 1, nk - 1):
        if (irout[k + 1] - irc) * (irout[k] - irc) <= 0:
            if abs(irout[k + 1] - irc)< 1E-9:
                x2 = k + 1
            elif  abs(irout[k] - irc)< 1E-9:
                x2 = k
            else:
                weight = abs(irout[k] - irc) / abs(irout[k + 1] - irout[k])
                x2 = k + weight
            break
   # print "x2", x2
    dzfwhm = vZsampling * (x2 -x1)
    irout *= maxirout


    return status, dzfwhm, irout










if '__main__' == __name__:
    impulse = zeros((2400))
   # status, dz, impulse =  NDAAC_ResolIr(1, array([1/5,1/5,1/5,1 /5.,1/5]), impulse)
 #   impulse = zeros((2400))
   # print NDAAC_ResolIr(1, 1/9*array([1,1,1,1,1,1,1,1,1]), impulse)

    impulse = zeros((2400))
   # print  NDAAC_ResolIr(1, array([-3,-2,-1,0,1,2,3]) / 28., impulse)
   # print NDAAC_ResolIr(1, ones((25)) / 25, impulse)
 #   status, dz, impulse = NDAAC_ResolIr(7.5, ones((25)) / 25, impulse)
    print("second")
    print(NDAAC_ResolIr(7.5, 1/2*array([-1,0,1]), impulse))
    impulse = zeros((2400))
    print(NDAAC_ResolIr(7.5 * 25, 1/2*array([-1,0,1]), impulse))
    






