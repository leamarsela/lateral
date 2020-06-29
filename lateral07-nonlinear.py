##################################
###Lateral Calculation for pile###
##################################

from math import pi
import numpy as np

#input
diameter = 0.5 #unit in m
modulus = 2.e7 #unit in kN/m2
inertia = pi * (diameter**4)/64. #unit in m4
Pt = 100. #unit in kN
constJ = 0.5
numsegmen = 100
numpoint = 16
stagept = 10

zo = [0., 3., 6., 11.]
zi = [3., 6., 11., 12.]
gamma = [10., 10., 10., 10.]
cu = [25., 25., 25., 25.]

numdata = len(zo)

numnodal = numsegmen + 1

numindexB = 2 * numnodal

deltaZ = zi[len(zi) - 1] / numsegmen


ndepth = []
for i in range(numnodal):
    if i == 0:
        ndepthI = 0
    else:
        ndepthI = deltaZ + ndepthI
    ndepth.append(round(ndepthI, 4))


ncu = []
ngamma = []
for i in range(numnodal):
    for j in range (len(cu)):
        if i == 0:
            ncui = cu[0]
            ngammai = gamma[0]
        elif i > 0 and i < numnodal - 1:
            if ndepth[i] > zo[j] and ndepth[i] < zi[j]:
                ncui = cu[j]
                ngammai = gamma[j]
        else:
            ncui = cu[len(cu) - 1]
            ngammai = gamma[len(gamma) - 1]
    ncu.append(ncui)
    ngamma.append(ngammai)


cui = []
gammai = []
for i in range(numnodal):
    if i == 0:
        templayer = ndepth[i+1] - ndepth[i]
        tempctimeslayer = ncu[i] * templayer
        tempgtimeslayer = ngamma[i] * templayer
        tempcui = tempctimeslayer / templayer
        tempgammai = tempgtimeslayer / templayer
    elif i > 0  and i < numnodal - 1:
        templayer = templayer + (ndepth[i+1] - ndepth[i])
        tempctimeslayer = tempctimeslayer + (ncu[i] * (ndepth[i+1] - ndepth[i]))
        tempgtimeslayer = tempgtimeslayer + (ngamma[i] * (ndepth[i+1] - ndepth[i]))
        tempcui = tempctimeslayer / templayer
        tempgammai = tempgtimeslayer / templayer
    else:
        tempcui = ncu[len(ncu) - 1]
        tempgammai = ngamma[len(ngamma) - 1]
    cui.append(round(tempcui, 2))
    gammai.append(round(tempgammai, 2))

def epsilon50(valcu):
    if valcu <= 48.:
        epsilon50 = 0.02
    elif valcu > 48.0 and valcu < 96.0:
        epsilon50 = 0.01
    else:
        epsilon50 = 0.005
    return epsilon50

def valy50(valepsilon50, valdiameter):
    return (2.5 * valepsilon50 * valdiameter)

def pult1(valgamma, valJ, valdepth, valcu, valdiameter):
    return ((3.0 + (valgamma * valdepth / valcu) + (valJ * valdepth / valdiameter)) * valcu * valdiameter)

def pult2(valcu, valdiameter):
    return (9.0 * valcu * valdiameter)

def valpu(valv, valpult, valy50):
    if valv/valy50 >= 8.0:
        pult = valpult
    else:
        if valv < 0.:
            pult = (0.5 * valpult * -abs(pow(valv/valy50, (1.0/3.0))))
        else:
            pult = (0.5 * valpult * (valv / valy50)**(1.0 / 3.0))
    return pult

def Dvalpu(valv, valpult, valy50):
    # return ((1 / 6.0) * (valpult / valy50) * (valv / valy50)**(-2.0 / 3.0))
    return ((1./6.) * (valpult / valy50) / ((pow(abs(valv)/valy50, (2.0/3.0)))))

def DDvalpu(valv, valpult, valy50):
    # return ((-1.0/9.0) * valpult / ((valy50**2) * (pow(abs(valv) / valy50, 5.0/3.0))))
    # return ((-valpult/9.0) / ((valy50**2) * (-abs(pow(valv/valy50, 5.0/3.0)))))
    return ((-valpult/9.0) / ((valy50**2) * ((pow(abs(valv)/valy50, 5.0/3.0)))))

def valko(valyo, valpu):
    return(valpu / valyo)


#method Newton Raphson
def PDPu(valv, valpult, valy50):
    return (valpu(valv, valpult, valy50) / Dvalpu(valv, valpult, valy50))


#method Raltson and Rabinowitz
def Ux(valv, valpult, valy50):
    return (valpu(valv, valpult, valy50) * Dvalpu(valv, valpult, valy50) / ((Dvalpu(valv, valpult, valy50)**2)-(valpu(valv, valpult, valy50) * DDvalpu(valv, valpult, valy50))))


#method Raltson and Rabinowitz multiple roots
def Uxi(dh, valv, valpult, valy50):
    valxmini = valv-dh
    valuxi = PDPu(valv, valpult, valy50)
    valuxmini = PDPu(valxmini, valpult, valy50)
    return (valuxi * (valxmini - valv) / (valuxmini - valuxi))



valpt = 0.
deltapt = Pt / stagept

outputvaly = []

for i in range(stagept):

    valpt = valpt + deltapt
    valC1 = 2.0 * valpt * (deltaZ)**3 / (modulus * inertia)


    kom = []
    for j in range(numnodal):
        try:
            # if i == 0:
            kom.append(Dvalpu(valy50(epsilon50(cui[j]), diameter)/1000, 
                              min(pult1(gammai[j], constJ, ndepth[j], cui[j], diameter), pult2(cui[j], diameter)), 
                              valy50(epsilon50(cui[j]), diameter)))

            # kom.append(valko(valy50(epsilon50(cui[j]), diameter)/100, 
            #                         valpu(valy50(epsilon50(cui[j]), diameter)/100, 
            #                             min(pult1(gammai[j], constJ, ndepth[j], cui[j], diameter), pult2(cui[j], diameter)),
            #                             valy50(epsilon50(cui[j]), diameter))))
            # # else:
            #     if (outputvaly[i-1][j] / valy50(epsilon50(cui[j]), diameter))  >= 8.0:
            #         kom.append(1.)
            #     else:
            # kom.append(Dvalpu((outputvaly[i-1][j]), 
            #                    min(pult1(gammai[j], constJ, ndepth[j], cui[j], diameter), pult2(cui[j], diameter)),
            #                    valy50(epsilon50(cui[j]), diameter)))
        except ZeroDivisionError:
            print('Zero Division')
    # print(kom)

    valA = []
    for j in range(numnodal):
        if j == 0:
            valA.append(0.)
        else:
            valA.append(kom[j] * (deltaZ**4) / (modulus * inertia))      
    valA.reverse()
    # print(valA)

    valB = []
    for k in range(numindexB):
        if k == 0:
            valB.append(2.0 / (valA[0] + 2.0))
        elif k == 1:
            valB.append(2 * valB[0])
        elif k == 2:
            valB.append(1.0 / (5.0 + valA[1] - 2.0 * valB[1]))
        else:
            if k % 2 != 0:
                m = (k-1)//2
                valB.append(valB[2*m] * (4.0 - valB[2 * m - 1]))
            elif k % 2 == 0:
                m = k // 2
                valB.append(1.0 / (6.0 + valA[m] - valB[2 * m -4] - valB[2 * m - 1] * (4.0 - valB[2 * m - 3])))
    # print(valB)

    valD = []
    for l in range(3):
        if l == 0:
            valD.append(1.0 / valB[2 * (numnodal - 1)])
        elif l == 1:
            valD.append(valD[0]*valB[2 * (numnodal - 1) + 1] - valB[2 * (numnodal -1) - 2] * (2.0 - valB[2 * (numnodal - 1) - 3]) - 2.0)
        elif l == 2:
            valD.append(valD[0] - valB[2 * (numnodal - 1) - 4] - valB[2 * (numnodal - 1) - 1] * (2.0 - valB[2 * (numnodal - 1) - 3]))
    # print(valD)

    valyi = []
    valyi.append(valC1 * (1 + valB[2 * (numnodal - 1) - 2]) / ((valD[2] * (1 + valB[2 * (numnodal - 1) - 2])) - (valD[1] * valB[2 * (numnodal - 1) - 1])))
    valyi.append(valB[2 * (numnodal - 1) - 1] * valyi[0] / (1 + valB[2 * (numnodal - 1) - 2]))
    valyi.append(valD[0] * valB[2 * (numnodal - 1) + 1] * valyi[1] - valD[0] * valyi[0])

    for j in range(numnodal - 2, -1, -1):
        valyi.insert(0, (-valB[2 * j] * valyi[1] + valB[2 * j + 1] * valyi[0]))

    valyi.insert(0, (2.0 * valyi[0] - valyi[1]))
    valyi.insert(0, (valyi[3] - 4.0 * valyi[2] + 4.0 * valyi[1]))    
    valyi.reverse()
    

    tol = 1.e-5
    iitermax = 100
    deltay = []
    finalvalv = 0.
    Dv = 0.
    for j in range(numnodal):
        
        valyij = valyi[j + 2]
        # print('initial', valyij)
        Pext = kom[j] * valyij
        # print('nilai ko', kom[j])
        iiter = 0
        error = 1.
        while error > tol:
            try:
                iiter += 1
                
                Dv = Dv + valyij
                # print(valyij)      
                # valnew = valyij - PDPu(valyij, 
                #                        min(pult1(gammai[j], constJ, ndepth[j], cui[j], diameter), pult2(cui[j], diameter)), 
                #                        valy50(epsilon50(cui[j]), diameter))

                # error = abs(valnew - valyij)
                # valyij = valnew
                # Pint = valpu(valyij + Dv, 
                #              min(pult1(gammai[j], constJ, ndepth[j], cui[j], diameter), pult2(cui[j], diameter)),
                #              valy50(epsilon50(cui[j]), diameter))

                Pint = valpu(Dv, 
                             min(pult1(gammai[j], constJ, ndepth[j], cui[j], diameter), pult2(cui[j], diameter)),
                             valy50(epsilon50(cui[j]), diameter))


                # if ((valyij + Dv) / valy50(epsilon50(cui[j]), diameter)) >= 8.0:
                #     Pext = Pint
                if (Dv / valy50(epsilon50(cui[j]), diameter)) >= 8.0:
                    Pext = Pint

                # dy = (Pext - Pint) / kom[j]

                dy = (Pext - Pint) / Dvalpu(Dv, 
                                            min(pult1(gammai[j], constJ, ndepth[j], cui[j], diameter), pult2(cui[j], diameter)),
                                            valy50(epsilon50(cui[j]), diameter))

                # valyij = valyij + dy
                Dv = Dv + dy

                error = abs(Pext - valpu(Dv, 
                                         min(pult1(gammai[j], constJ, ndepth[j], cui[j], diameter), pult2(cui[j], diameter)),
                                         valy50(epsilon50(cui[j]), diameter)))
                # print(error, Pext, kom[j], Dv, iiter, valyij)
                # print('nilai Pext', Pext)
                # print('nilai Po', valpu(valyij + Dv, 
                #                          min(pult1(gammai[j], constJ, ndepth[j], cui[j], diameter), pult2(cui[j], diameter)),
                #                          valy50(epsilon50(cui[j]), diameter)))
                # print('error', error)

                if iiter == iitermax:
                    break
                
            except ZeroDivisionError:
                print('error')
        
        finalvalv = finalvalv + Dv
        valyij = 0.
        # print('nilai dv', Dv)
        # valyij = valyij + Dv
        # print('initial valyij', valyij, kom[j])
        # print('nilai vlyij', valyij, Dv)
        # Dv = 0.
        deltay.append(finalvalv)
    outputvaly.append(deltay)

print(outputvaly)



##untuk cetak hasil
# from openpyxl import Workbook

# wb = Workbook()
# ws = wb.active

# valpt = 0
# for i in range(stagept):
#     valpt = valpt + deltapt
#     ws.cell(2 + i, 1).value = valpt

# for i in range(len(outputvaly)):
#     for j in range(numnodal):
#         ws.cell(2 + i, 2 + j).value = outputvaly[i][j]

# namefile = '10pileversi' + str(numsegmen) + str(zi[len(zi) - 1]) + str(cu[0]) + str(diameter*10) + '.' + 'xlsx'

# wb.save('./results/' + namefile)


    