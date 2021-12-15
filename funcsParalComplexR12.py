import matplotlib.pyplot as plt
from datetime import datetime
from multiprocessing import Pool
import scipy.special as sspecial
import mpmath
from mpmath import mp
import numpy as np
mp.dps=30

#this script computes potential 4 V=x^{2}-igx^{5}, region I-II

def ret5AdjacentPairs(g,E):
    """

    :param g: const
    :param E: trial eigenvalue
    :return: 5 adjacent pairs of roots, the first root x2 has smaller angle than the second root x1,
    return [[x2, x1]]
    """
    coefs = [-1j * g, 0, 0, 1, 0, -E]
    rootsAll=mpmath.polyroots(coefs)

    rootsSortedByAngle =sorted(rootsAll,key=mpmath.arg)
    rst=[]
    length=len(rootsSortedByAngle)
    for j in range(0,length):
        rst.append([rootsSortedByAngle[j],rootsSortedByAngle[(j+1)%length]])

    return rst


def ret5SeparatedPairs(g,E):
    """

    :param g: const
    :param E: trial eigenfunction
    :return: 5 pairs of roots separated by another root, the first root x2 has smaller angle than the
    second root x1, return [[x2, x1]]
    """
    coefs = [-1j * g, 0, 0, 1, 0, -E]
    rootsAll = mpmath.polyroots(coefs)

    rootsSortedByAngle = sorted(rootsAll, key=mpmath.arg)
    rst = []
    length = len(rootsSortedByAngle)
    for j in range(0, length):
        rst.append([rootsSortedByAngle[j], rootsSortedByAngle[(j + 2) % length]])

    return rst


def f(z, g, E):
    '''
    :param g: const
    :param z: point on x2x1
    :param E: trial eigenvalue
    :return: f value
    '''
    return mp.sqrt(1j * g * z ** 5 - z ** 2 + E)


def fBranchOther(z,g,E):
    '''

    :param z: point on x2x1
    :param g: const
    :param E: trial eigenvalue
    :return: f value on another branch
    '''
    return -mp.sqrt(1j * g * z ** 5 - z ** 2 + E)

def integralQuadrature(g,E,x1,x2):
    '''

    :param g: const
    :param E: trial eigenvalue
    :param x1: ending point
    :param x2: starting point
    :return:
    '''
    a1 = mpmath.re(x1)
    b1 = mpmath.im(x1)

    a2 = mpmath.re(x2)
    b2 = mpmath.im(x2)
    slope = (b1 - b2) / (a1 - a2)
    gFunc = lambda y: f(y + 1j * (slope * (y - a2) + b2), g, E)
    return (1 + 1j * slope) * mpmath.quad(gFunc, [a2, a1])



def integralQuadratureAnotherBranch(g,E,x1,x2):
    '''

    :param g: const
    :param E: trial eigenvalue
    :param x1: ending point
    :param x2: starting point
    :return:
    '''
    a1 = mpmath.re(x1)
    b1 = mpmath.im(x1)

    a2 = mpmath.re(x2)
    b2 = mpmath.im(x2)
    slope = (b1 - b2) / (a1 - a2)
    gFunc = lambda y: fBranchOther(y + 1j * (slope * (y - a2) + b2), g, E)
    return (1 + 1j * slope) * mpmath.quad(gFunc, [a2, a1])



def eqnFiveAdjacentPairs(EIn,n,g):
    """
    computes adjacent pairs
    :param EIn:
    :param n:
    :param g:
    :return:
    """
    E=EIn
    adjPairsAll = ret5AdjacentPairs(g, E)
    retValsCis = []  # in the order x2, x1
    retValsTrans = []  # in the order x1,x2
    retValsCisAnother = []  # in the order x2, x1, another branch
    retValsTransAnother = []  # in the order x1, x2, another branch

    ###############################
    # fill cis

    # for pairTmp in adjPairsAll:
    #     x2Tmp, x1Tmp = pairTmp
    #     intValTmp = integralQuadrature(g, E, x1Tmp, x2Tmp)
    #     rstTmp = intValTmp - (n + 1 / 2) * mp.pi
    #     retValsCis.append(rstTmp)
    # # fill trans
    # for pairTmp in adjPairsAll:
    #     x2Tmp, x1Tmp = pairTmp
    #     intValTmp = integralQuadrature(g, E, x2Tmp, x1Tmp)
    #     rstTmp = intValTmp - (n + 1 / 2) * mp.pi
    #     retValsTrans.append(rstTmp)
    ########################################
    ########################################
    # # #fill cis another
    for pairTmp in adjPairsAll:
        x2Tmp, x1Tmp = pairTmp
        intValTmp = integralQuadratureAnotherBranch(g, E, x1Tmp, x2Tmp)
        rstTmp = intValTmp - (n + 1 / 2) * mp.pi
        retValsCisAnother.append(rstTmp)

    # #fill trans another
    for pairTmp in adjPairsAll:
        x2Tmp, x1Tmp = pairTmp
        intValTmp = integralQuadratureAnotherBranch(g, E, x2Tmp, x1Tmp)
        rstTmp = intValTmp - (n + 1 / 2) * mp.pi
        retValsTransAnother.append(rstTmp)
    #################################################
    retCombined = retValsCis + retValsTrans + retValsCisAnother + retValsTransAnother
    retSorted = sorted(retCombined, key=mpmath.fabs)
    root0 = retSorted[0]
    return root0



def computeOneSolutionWith5AdjacentPairs(inData):
    """

        :param inData: [n, g, Eest]
        :return: [n, g, E]
        """
    n, g, Eest = inData
    Eest*=(1+0j)
    try:
        E=mpmath.findroot(lambda EVal:eqnFiveAdjacentPairs(EVal,n,g),Eest,solver="muller",maxsteps=100,tol=1e-10)
        return [n, g, E]
    except ValueError as e:
        # print("error message:"+str(e))
        return []


def vecEqnFiveAdjacentPairs(ERe,EIm,n,g):
    """

    :param ERe:
    :param EIm:
    :param n:
    :param g:
    :return:
    """
    E=ERe+1j*EIm
    adjPairsAll = ret5AdjacentPairs(g, E)
    retValsCis = []  # in the order x2, x1
    retValsTrans = []  # in the order x1,x2
    retValsCisAnother = []  # in the order x2, x1, another branch
    retValsTransAnother = []  # in the order x1, x2, another branch

    ###############################
    # fill cis

    for pairTmp in adjPairsAll:
        x2Tmp, x1Tmp = pairTmp
        intValTmp = integralQuadrature(g, E, x1Tmp, x2Tmp)
        rstTmp = intValTmp - (n + 1 / 2) * mp.pi
        retValsCis.append(rstTmp)
    # fill trans
    for pairTmp in adjPairsAll:
        x2Tmp, x1Tmp = pairTmp
        intValTmp = integralQuadrature(g, E, x2Tmp, x1Tmp)
        rstTmp = intValTmp - (n + 1 / 2) * mp.pi
        retValsTrans.append(rstTmp)
    ########################################
    ########################################
    # # #fill cis another
    # for pairTmp in adjPairsAll:
    #     x2Tmp, x1Tmp = pairTmp
    #     intValTmp = integralQuadratureAnotherBranch(g, E, x1Tmp, x2Tmp)
    #     rstTmp = intValTmp - (n + 1 / 2) * mp.pi
    #     retValsCisAnother.append(rstTmp)
    #
    # # #fill trans another
    # for pairTmp in adjPairsAll:
    #     x2Tmp, x1Tmp = pairTmp
    #     intValTmp = integralQuadratureAnotherBranch(g, E, x2Tmp, x1Tmp)
    #     rstTmp = intValTmp - (n + 1 / 2) * mp.pi
    #     retValsTransAnother.append(rstTmp)
    #################################################
    retCombined = retValsCis + retValsTrans + retValsCisAnother + retValsTransAnother
    retSorted = sorted(retCombined, key=mpmath.fabs)
    root0 = retSorted[0]
    return np.real(root0),np.imag(root0)


def vecComputeOneSolutionWith5AdjacentPairs(inData):
    """

    :param inData: [n, g, Eest]
    :return: [n,g,E]
    """
    n, g, Eest = inData
    EReEst=np.real(Eest)
    EImEst=np.imag(Eest)
    func=lambda ERe, EIm: vecEqnFiveAdjacentPairs(ERe,EIm,n,g)
    try:
        EReVal,EImVal=mpmath.findroot(func,(EReEst,EImEst),solver="halley",maxsteps=100,tol=1e-10)

        E=EReVal+1j*EImVal
        return [n,g,E]
    except ValueError as e:
        # print("error message:"+str(e))
        return []
