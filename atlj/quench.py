#!/usr/bin/env python3
from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
from builtins import object
import numpy as np
import pickle
from pele.optimize import LBFGS_CPP
from pele.optimize import MixedOptimizer
import itertools
import os

# def lbfgs_cpp_steprecorder(coords, pot, niter, **kwargs):
#     lbfgs = LBFGS_CPP(coords, pot, **kwargs)
#     # also returnp
#     elist, gradlist = [], []
#     coordslist = []
#     for i in range(int(niter)):
#         # if lbfgs.stop_criterion_satisfied():
#         #     break
#         res = lbfgs.one_iteration()
#         coords = res.coords
#         e, grad = pot.getEnergyGradient(coords)
#         elist.append(e)
#         gradlist.append(grad)
#         coordslist.append(coords)
#     history = [np.array(elist), np.array(gradlist), np.array(coordslist)]
#     return copy.deepcopy(res), history


class PlaneTransform(object):
    def __init__(self, alpha):
        self.alpha = alpha
        self.vp = np.array([1., 1., 1.]) #the vector perpendicular to the plane
        v1 = np.array([0.,0,1]) # an arbitrary vector
        #n1 and n2 are orthogonal unit vectors which define the plane
        self.n1 = np.cross(self.vp,v1)
        self.n2 = np.cross(self.vp, self.n1)
        #normalize them
        self.n1 /= np.linalg.norm(self.n1)
        self.n2 /= np.linalg.norm(self.n2)
        #self.vp /= np.linalg.norm(self.vp)
        #flip the axes, to get them to be like david's
        self.n1, self.n2 = self.n2, self.n1
        print(self.n1)
        print(self.n2)
        print(np.cross(self.n1,self.n2))
        
    def projectOntoPlane(self,Xin):
        #find components parallel to n1 and n2
        #X = Xin-self.alpha
        Xp = np.zeros(2)
        Xp[0] = np.dot(X,self.n1)
        Xp[1] = np.dot(X,self.n2)
        return Xp
    
    def PlaneToXYZ(self, Xp):
        return Xp[0]*self.n1 + Xp[1]*self.n2 + self.alpha*self.vp

def xyz2internal(X):
    """
    Cartesian to internal coords.
    """

    X = np.reshape(X, [3*3])

    R = np.zeros([3,3])

    R[0] = X[0:3]-X[3:6]
    R[1] = X[3:6]-X[6:9]
    R[2] = X[0:3]-X[6:9]

    return R

def internal2xyz(R):
    """
    assert the three particles lie in the plane z=0, that X2 is at the origin,
    and that X1 is along the x axis
    """
    X = np.zeros([3,3])
    #X2 = np.zeros(3)
    X[0,:] = np.array([R[0], 0.,0.])
    #from the cosine rule
    costheta123 = old_div(-(R[2]**2 - R[0]**2 - R[1]**2), (2*R[1]*R[0]))
    if costheta123 > 1:
        print("warning: costheta123 > 1:", costheta123)

    X3x = R[1] * costheta123
    X3y = R[1] * np.sqrt(1 - costheta123**2)
    X[2,:] = np.array([X3x, X3y, 0.])
    
    if False: #do some testing
        i = 0; j = 1; k = 0
        r = np.linalg.norm(X[i,:]-X[j,:])
        if r - R[k] > 1e-4:
            print("error: internal2xyz", k)
            print(X[i,:], X[j,:], r, R[k], np.abs(r-R[k]), costheta123, R)
            print(X)
            exit(1)

        i = 1; j = 2; k = 1
        r = np.linalg.norm(X[i,:]-X[j,:])
        if r - R[k] > 1e-4:
            print("error: internal2xyz", k)
            print(X[i,:], X[j,:], r, R[k], np.abs(r-R[k]), costheta123, R)
            print(X)
            exit(1)

        i = 2; j = 0; k = 2
        r = np.linalg.norm(X[i,:]-X[j,:])
        if r - R[k] > 1e-4:
            print("error: internal2xyz", k)
            print(X[i,:], X[j,:], r, R[k], np.abs(r-R[k]), costheta123, R)
            print(X)
            exit(1)

    return X


def getCoords(alpha = 2.**(1./6), nmesh = 200 ):
    #from pele.potentials.ATLJ import ATLJ
    #atlj = ATLJ(Z = Z)

    #take a cut along the plain perpendicular to r12 = r13 = r23
    #i.e. r12 + r13 + r23 = const = alpha
    #alpha = 2.#2.**(1./6)
    
    #let xp, yp be the axes in this plane
        
    xmax = alpha
    xmin = -alpha
    ymax = alpha
    ymin = -alpha
    
    #ymin = alpha/7.
    #ymax = 3*alpha/7
    #xmin = 2*alpha/7
    #xmax = 4*alpha/7
    
    nmesh = nmesh
    dx = old_div((xmax-xmin),nmesh)
    dy = old_div((ymax-ymin),nmesh)
    xpmesh = np.array([xmin + j*dx for j in range(nmesh)])
    ypmesh = np.array([ymin + j*dy for j in range(nmesh)])
        
    tform = PlaneTransform(alpha)
    
    if False:
        for i in range(100):
            x = np.random.uniform(0,2,3)
            xp = tform.projectOntoPlane(x)
            xnew = tform.PlaneToXYZ(xp)
            print(x[0], x[1], x[2], xp[0], xp[1], xnew[0], xnew[1], xnew[2])

    # from pele.printing.print_atoms_xyz import printAtomsXYZ as printxyz
    coordslist = []
    rej_triangle = 0
    rej_close = 0
    for xp in xpmesh:
        for yp in ypmesh:
            Xp = np.array([xp,yp])
            
            Rinternal = tform.PlaneToXYZ(Xp)
            
            #check if it is physical
            reject = False
            for i in range(3):
                if Rinternal[i] + Rinternal[i-1] < Rinternal[i-2]:
                    reject = True
                    rej_triangle += 1
                    break
            if (Rinternal < 0.9).any():
                rej_close += 1 
                reject=True
            if reject:
                #print "rejecting", xp, yp, Rinternal
                continue
            
            coords = internal2xyz(Rinternal)
            coords = np.reshape(coords, [3*3])
            coordslist.append( (coords, np.array([xp,yp])))
    
    print(rej_triangle, "configurations rejected because of the triangle inequality")
    print(rej_close, "configurations rejected because the atoms were too close")
    # with open("starting.xyz","w") as fout:
    #     # for coords,Xp in coordslist:
    #     #     printxyz(fout, coords)

    return coordslist






# from pele.optimize._quench import lbfgs_cpp as quench
from pele.optimize._quench import steepest_descent
from pele.optimize import GradientDescent_CPP
# from pele.optimize._quench import cg as quench
#from pele.optimize.quench import quench as quench
#from pele.optimize.quench import bfgs as quench
# from pele.optimize._quench import lbfgs_py as quench
#from pele.optimize.quench import mylbfgs as quench
from pele.optimize import ModifiedFireCPP
from pele.potentials import ATLJCPP as ATLJ
# from pele.optimize import LBFGS_MPFR_CPP
from pele.optimize._quench import ode_scipy_naive
from pele.optimize._quench_jl import ode_julia_naive






def testCoords(coordslist, Z=2., opt_class=LBFGS_CPP):
    
    #quench.gtol = 1.e-5

    pot = ATLJ(Z = Z)

    finish = []
    count = 0
    ntot = len(coordslist)
    optimizer_initialized = False
    for coords, Xp in coordslist:
        if optimizer_initialized is False:
            optimizer = opt_class(coords, pot, T=5)
            optimizer_initialized = True
        else:
            optimizer.reset(coords)
        if count % 1000 == 0:
            print("quenching", count, "out of", ntot)
        count += 1
        initen = pot.getEnergy(coords)
        ret = optimizer.run(5)
        coordsnew = ret.coords
        print(ret.nfev)
        finish.append((Xp, OrderParam(xyz2internal(coordsnew)), initen))
    with open('finish_radford_main.pkl', 'wb') as f:
        pickle.dump(finish, f)
    print("this works")
    # from pele.utils.xyz import printAtomsXYZ as printxyz
    # with open("finishradford.xyz","w") as fout:
    #    for ret in finish:
    #        coords= ret[2]
    #        e = ret[3]
    #        printxyz(fout, coords, line2=str(e))

    return finish



# def testCoordscompare(coordslist, Z = 2., opt_class=LBFGS_MPFR_CPP, precision=2):
#     """ Compares to the steepest descent version of the potential
#     """
#     pot = ATLJ(Z = Z)

#     finish = []
#     count = 0
#     ntot = len(coordslist)
#     optimizer_initialized = False
#     for coords, Xp in coordslist:
#         if optimizer_initialized is False:
#             optimizer = opt_class(coords, pot, prec=precision)
#             optimizer_initialized = True
#         else:
#             optimizer.reset(coords)
#         if count % 1000 == 0:
#             print("quenching", count, "out of", ntot)
#         initen = pot.getEnergy(coords)
#         ret = optimizer.run(2000)
#         coordsnew = ret.coords
#         print(ret.nfev)
#         finish.append((Xp, OrderParam(xyz2internal(coordsnew)), initen))
#         count += 1

#     with open('finish_radford_main.pkl', 'wb') as f:
#         pickle.dump(finish, f)
#     print("this works")
#     # from pele.utils.xyz import printAtomsXYZ as printxyz
#     # with open("finishradford.xyz","w") as fout:
#     #    for ret in finish:
#     #        coords= ret[2]
#     #        e = ret[3]
#     #        printxyz(fout, coords, line2=str(e))

#     return finish


def testCoordsquench(coordslist, Z = 2., quench=ode_julia_naive, precision=10, opt_class=LBFGS_CPP):
    """ Compares to the steepest descent version of the potential
    """
    pot = ATLJ(Z = Z)
    finish = []
    count = 0
    ntot = len(coordslist)
    optimizer_initialized = False
    correctminima = []
    orderparamlist = []
    nfevlist = []
    nstepslist = []
    for coords, Xp in coordslist:
        if count % 1 == 0:
            print("quenching", count, "out of", ntot)
        count += 1
        initen = pot.getEnergy(coords)
        ret = quench(coords, pot)
        coordsnew = ret.coords
        nfevlist.append(ret.nfev)
        nstepslist.append(ret.nsteps)
        print(ret.nfev, 'function evaluations')
        op1 = OrderParam(xyz2internal(coordsnew))
        orderparamlist.append(op1)
        finish.append((Xp, OrderParam(xyz2internal(coordsnew)), initen))
    print(np.average(nfevlist))
    print(np.average(nstepslist))
    np.savetxt('/home/praharsh/Dropbox/research/bv-libraries/orderparamliststeepest.txt', orderparamlist, delimiter=',')
    # percent = np.average(np.array(correctminima, dtype=float))
    # print(percent)
    print("this works")
    # from pele.utils.xyz import printAtomsXYZ as printxyz
    # with open("finishradford.xyz","w") as fout:
    #    for ret in finish:
    #        coords= ret[2]
    #        e = ret[3]
    #        printxyz(fout, coords, line2=str(e))
    return finish



def testCoordsquenchl(coordslist, Z = 2., quench=steepest_descent, precision=23, opt_class=MixedOptimizer, T=1, M=1):
    """ Tests 
    """
    pot = ATLJ(Z=Z)
    finish = []
    count = 0
    ntot = len(coordslist)
    optimizer_initialized = False
    correctminima = []
    opsteep = np.loadtxt('/home/praharsh/Dropbox/research/bv-libraries/orderparamliststeepest.txt', delimiter=',')
    nfevlist = []
    niterlist = []
    oplist = []
    # [16690:] not working nmesh=700 backtracking
    for coords, Xp in coordslist:
        if count % 1 == 0:
            print("quenching", count, "out of", ntot)
        # optimizer = opt_class(coords, pot, M=M, T=T, iprint=-1, maxstep=1, tol=1e-4)
        # for steepest descent
        optimizer = opt_class(coords, pot, tol=1e-4, T=10)
        # optimizer_initialized = True
        # optimizer = opt_class(coords, pot, M=M, T=T, iprint=-1, maxstep=1, tol = 1e-4)
        # optimizer_initialized = True
        initen = pot.getEnergy(coords)
        optimizer.set_H0(1e-6)
        ret2 = optimizer.run(100)
        op1 = opsteep[count]
        # op1 = 0
        coordsnew2 = ret2.coords
        op2 = OrderParam(xyz2internal(coordsnew2))
        print(coordsnew2)
        correctminima.append(op1 == op2)
        # print(op1 == op2)
        oplist.append(op2)
        finish.append((Xp, OrderParam(xyz2internal(coordsnew2)), initen))
        nfevlist.append(ret2.nfev)
        print(ret2.nfev)
        niterlist.append(ret2.nsteps)
        count += 1

    th = 1-np.average(np.array(correctminima, dtype=float))
    np.savetxt('/home/praharsh/Dropbox/research/bv-libraries/nfevlistM4precon.txt', nfevlist, delimiter=',')
    np.savetxt('/home/praharsh/Dropbox/research/bv-libraries/orderparamliststeepesdescentnew.txt', oplist, delimiter=',')
    print(th, "percentage of wrong minima")
    # from pele.utils.xyz import printAtomsXYZ as printxyz
    # with open("finishradford.xyz","w") as fout:
    #    for ret in finish:
    #        coords= ret[2]
    #        e = ret[3]
    #        printxyz(fout, coords, line2=str(e))
    print(np.average(niterlist), "number of iterations")
    print(np.average(nfevlist), "number of function evaluations")
    return finish, np.average(nfevlist), th





def quenchsingle(coordsXp, pot, opt_class):
    coords, Xp = coordsXp
    print(Xp)
    initen = pot.getEnergy(coords)
    optimizer = opt_class(coords, pot, M=1, T=2, iprint=-1, maxstep=0.1)
    ret = optimizer.run(30)
    coordsnew = ret.coords
    nfevel = ret.nfev
    return ((Xp, OrderParam(xyz2internal(coordsnew)), initen), nfevel)



def OrderParam(R):
    
    costheta1 = np.abs(np.dot(old_div(R[0],np.linalg.norm(R[0])),old_div(R[1],np.linalg.norm(R[1]))))
    costheta2 = np.abs(np.dot(old_div(R[1],np.linalg.norm(R[1])),old_div(R[2],np.linalg.norm(R[2]))))
    costheta3 = np.abs(np.dot(old_div(R[0],np.linalg.norm(R[0])),old_div(R[2],np.linalg.norm(R[2]))))
    
    ##### Failed quench, not linear or triangular or too far apart.
    struct = 5.0 
    
    ##### Check distances betweet particles.
    if np.abs(np.linalg.norm(R[0]))+np.abs(np.linalg.norm(R[1]))+np.abs(np.linalg.norm(R[2])) < 4.4*2.**(1./6.):
        
        ##### It's a TRIANGLE.
        if 0.4 < costheta1 < 0.6 and 0.4 < costheta2 < 0.6 and 0.4 < costheta3 < 0.6: 
            struct = 1.0 
            
        ##### It's LINEAR.
        elif 0.9 < costheta1 and 0.9 < costheta2 and 0.9 < costheta3: 
            
            if np.abs(np.linalg.norm(R[0])-np.linalg.norm(R[2])) < 0.1:
                struct = 2.0 
                
            elif np.abs(np.linalg.norm(R[1])-np.linalg.norm(R[2])) < 0.1:
                struct = 3.0 
                
            else:
                struct = 4.0
    return struct

if __name__ == "__main__":
    picklef = "finish.pkl"
    Z = 2
    Re = 2.**(1./6.)
    alpha = np.sqrt(3.)*Re
    nmesh = 50
    
    coordslist = getCoords(alpha = alpha, nmesh=nmesh)
    pot = ATLJ(Z=Z)
    path = '/home/praharsh/Dropbox/research/bv-libraries/tests/data'

    # M = 4
    # T = 1
    
    # try:
    #     with open(picklef, "r") as fin:
    #         finish = pickle.load(fin)
    #     print("loading data from pickle file", picklef)
    # except:
    coordslist = getCoords(alpha=alpha, nmesh=nmesh)
    # finish, nfevlist, correctminima = testCoordsquenchl(coordslist, Z=Z, T=T, M=M)
    finish = testCoordsquench(coordslist, Z=Z)
    # with open(picklef, "w") as fout:
    # pickle.dump(finish, fout)
    import sys
    
    # ##########################
    # #####    PLOTTING    #####
    # ##########################

    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib import colors
    from matplotlib import rc
    
    ##### Use LaTeX fonts.
    rc('text', usetex=True) 
    rc('font',**{'family':'serif','serif':['Computer Modern']})
        
    xmax = alpha
    xmin = -alpha
    ymax = alpha
    ymin = -alpha
    
    #ymin = alpha/7
    #ymax = 3*alpha/7
    #xmin = 2*alpha/7
    #xmax = 4*alpha/7        
    xfact = old_div(nmesh,(xmax-xmin))
    yfact = old_div(nmesh,(ymax-ymin))
    
    grid = np.zeros( [nmesh, nmesh] )

    for data in finish:
        x = data[0][0]
        xin = int((x-xmin)*xfact-0.5)
        y = data[0][1]
        yin = int((y-ymin)*yfact-0.5)
        z = data[1] # 1 is order param, 2 is init energy.
        grid[xin,yin] = z
    
    ##### Define colors.
    cmap = colors.ListedColormap([(1.0,1.0,1.0), (52./255.,23./255.,0.0), (42./255.,95./255.,153./255.), (1.0,127./255.,0.0), (1.0,1.0,0.0),(127./255.,0.0,127./255.)])
    bounds=[-0.5,0.5,1.5,2.5,3.5,4.5,5.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    plt.imshow(grid,interpolation='nearest',cmap=cmap,norm=norm)
    #plt.imshow(grid)
    #plt.colorbar()
    #plt.gca().axison = False
    plt.xticks([])
    plt.yticks([])
    #plt.grid(True)
    #plt.show()
    filedir = "morethuente"
    os.makedirs(filedir, exist_ok=True)
    # plt.savefig( filedir + '/' + str(M) + 'T' + str(T) + 'nfev_' + str(nfevlist) + 'cm_' + str(correctminima) +'.pdf', bbox_inches='tight')
    plt.savefig('out.pdf', bbox_inches='tight')
