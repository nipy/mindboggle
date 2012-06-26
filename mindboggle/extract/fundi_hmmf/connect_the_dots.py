#!/usr/bin/python
"""
Connect the dots.

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com  (original Matlab code)
Arno Klein  .  arno@mindboggle.info  (translated to Python)

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

from simple_point_test import simple_point_test

#--------------------
# Compute probability
#--------------------
def p = compute_probability(wl,li,q,qneighs,wneighs):
    """
    Compute probability
    """
    p = -1 * (np.sqrt((wl-li)^2)*q + (wneighs * sum((q - qneighs).^2)))

#-----------------------
# Test for simple points
#-----------------------
def simple_point_test(faces, ind, values):
    import numpy as np
    from find_neighbors import find_neigbhors

    thr = 0.5
    neighs = find_neighbors(faces,ind)
    neighvals = values(neighs)
    nOutside = sum(neighvals <= thr)
    nInside = sum(neighvals > thr)

    if (nOutside == 0 || nInside == 0)
        sp = 0
        return

    if (nOutside == 1)
        sp = 1
        return

    if (nInside == 1)
        sp = 1
        return

    neighsInside = neighs(neighvals > thr)
    nSets = size(neighsInside,1)

    sets = zeros(nSets,20)
    setLabels = np.transpose(range(1,nSets+1))

    values(ind) = -1

    for i in range(nSets):
        currNeighs = find_neighbors(faces,neighsInside(i))
        currNeighs = currNeighs(values(currNeighs)>thr)
        sets(i,1:size(currNeighs,1)) = np.transpose(currNeighs)
        sets(i,size(currNeighs,1) + 1) = neighsInside(i)

    change = 1
    while (change > 0):
        change = 0
        for i in range(nSets):
            for j in range(nSets):
               if (j>i && setLabels(i) != setLabels(j)):
                    iInds = np.unique(sets(i,:))
                    iInds = iInds(iInds > 0)
                    iSize = size(iInds,2)

                    jInds = np.unique(sets(j,:))
                    jInds = jInds(jInds > 0)
                    jSize = size(jInds,2)

                    totInds = unique([iInds jInds])
                    totSize = size(totInds,2)

                    if (totSize < iSize + jSize):

                        minLabel = min(setLabels(j),setLabels(i))

                        setLabels(i) = minLabel
                        setLabels(j) = minLabel

                        change = 1

    numberOfSeparateSets = size(unique(setLabels),1)
    if (numberOfSeparateSets < 2):
        sp = 1
    else:
        sp = 0

    return sp

#=================
# Connect the dots
#=================
def Q = connectTheDotsV3(L, initValues, vertices, faces, dots, nList, shortInds):
    """
    Connect the dots.
    """
    likelihoodLimit = .5
    stepSize = .05
    wneighs = .4
    wl = 1.1
    multModInit = .02
    multMod = multModInit

    Q = zeros(length(L),1)
    Q(initValues>likelihoodLimit) = initValues(initValues>likelihoodLimit)

    Q(dots > .5) = 1

    if (sum(Q>0) < 2):
        return

    m = length(L)

    print('Initial candidates:')
    print(sum(Q>0))

    iterationCount = 0

    QNEW = Q

    pPrev = zeros(size(Q))

    for i = 1:size(pPrev,1):
        neighs = nList(shortInds == i,:)
        neighs = neighs(neighs > 0)
        pPrev(i) = compute_probability(wl,L(i),Q(i),Q(neighs),wneighs)

    totalPr = sum(pPrev(L>0))
    totalPr2 = totalPr - 10
    prevFundPoints = Inf

    % endFlag is used to not stop the iteration at the first occurance of no
    % change
    endFlag = 3

    while (endFlag > 0 && iterationCount < 350):

        iterationCount = iterationCount + 1
        downCounter = 0
        upCounter = 0

        multMod = multModInit + iterationCount*.001

        for i= 1:m:
            if (L(i) > 0 && Q(i) > .01):

                q = Q(i)

                qDown = max(q - stepSize, 0)

                neighs = nList(shortInds == i,:)
                neighs = neighs(neighs > 0)

                gPDown = compute_probability(wl,L(i),qDown,Q(neighs),wneighs)

                grDown = (gPDown - pPrev(i))/stepSize
                grDown = grDown * multMod

                if (grDown > 0):

                    if (q - grDown <= likelihoodLimit && q > likelihoodLimit):
                        if (dots(i) > .5):
                            c = 0
                        else:
                            QConnectivity = QNEW
                            QConnectivity(i) = q - grDown
                            c = simple_point_test(faces, i, QConnectivity)
                    else:
                        c =1
                    if (c == 1):
                        QNEW(i) = max(q - grDown, 0)
                        pPrev(i) = compute_probability(wl,L(i),QNEW(i),Q(neighs),wneighs)
                        downCounter = downCounter + 1
                else:
                    if (q - grDown > likelihoodLimit && q <= likelihoodLimit):
                        QConnectivity = QNEW
                        QConnectivity(i) = q - grDown
                        QConnectivity = 1+(QConnectivity *-1)
                        c = simple_point_test(faces, i, QConnectivity)
                    else:
                        c = 1
                    if (c == 1):
                        QNEW(i) = min(q - grDown, 1)
                        pPrev(i) = compute_probability(wl,L(i),QNEW(i),Q(neighs),wneighs)
                        upCounter = upCounter + 1

        diff = np.sqrt(sum((Q-QNEW).^2))

        totalPr = sum(pPrev(L>0))
        totalPrDiffPerPoint = (totalPr - totalPr2)/sum(L>0)

        currFundPoints = sum(Q>0.5)

        if (totalPrDiffPerPoint < 0.0001 && (currFundPoints - prevFundPoints) == 0):
            endFlag = endFlag - 1

        totalPr2 = totalPr
        print([diff totalPrDiffPerPoint])
        prevFundPoints = currFundPoints

        Q = QNEW
        print(['Iteration: ' num2str(iterationCount) ' ** Up and down movement: ' num2str(upCounter) ', ' num2str(downCounter) ' ** Current fundus points (q > .5): ' num2str(currFundPoints)])
