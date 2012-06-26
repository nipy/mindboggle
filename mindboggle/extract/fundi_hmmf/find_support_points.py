function P = findSupportPoints(fv,values,Umin)

thr = 5;

P = zeros(size(values));
checkList = P;

[maxVal maxInd] = max(values);
checkList(maxInd) = 1;
P(maxInd) = 1;
values(maxInd) = -1;

while (min(checkList) == 0 && maxVal > .5)
    
    [maxVal maxInd] = max(values);
    checkList(maxInd) = 1;
    values(maxInd) = -1;
    
    currPVert = fv.vertices(P>0,:);
    currU = Umin(:,P>0);    
    i = 0;
    found = 0;
    while(i < size(currPVert,1) && found ==0)
        i = i + 1;        
        D = compEucDist(currPVert(i,:),fv.vertices(maxInd,:));
        
        if (D >= thr && D<8)
            D = compDirectionalDist(currPVert(i,:),fv.vertices(maxInd,:),currU(:,i));
        end
        if (D < thr)
            found = 1;
        end
    end
    if (found == 0)        
        P(maxInd) = 1;
    end
end
