function Q = connectTheDotsV3(L,initValues,fv,dots,nList,shortInds)

%%%%%%%%%
likelihoodLimit = .5;
stepSize = .05;
%%%%%%%%%
wneighs = .4;
wl =1.1;
multModInit = .02;
multMod = multModInit;
%%%%%%%%%

Q = zeros(length(L),1);
Q(initValues>likelihoodLimit) = initValues(initValues>likelihoodLimit);

Q(dots > .5) = 1;

if (sum(Q>0) < 2)
    return;
end

m = length(L);

disp('Initial candidates:');
disp(sum(Q>0));

iterationCount = 0;

QNEW = Q;

pPrev = zeros(size(Q));

for i = 1:size(pPrev,1)
    neighs = nList(shortInds == i,:);
    neighs = neighs(neighs > 0);
    pPrev(i) = compProb(wl,L(i),Q(i),Q(neighs),wneighs);
end

totalPr = sum(pPrev(L>0));
totalPr2 = totalPr - 10;
prevFundPoints = Inf;

% endFlag is used to not stop the iteration at the first occurance of no
% change
endFlag = 3;

while (endFlag > 0 && iterationCount < 350)
    
    iterationCount = iterationCount + 1;
    downCounter = 0;
    upCounter = 0;
        
    multMod = multModInit + iterationCount*.001;
    
    for i= 1:m
        if (L(i) > 0 && Q(i) > .01)
            
            q = Q(i);
                                               
            qDown = max(q - stepSize, 0);            
            
            neighs = nList(shortInds == i,:);
            neighs = neighs(neighs > 0);
            
            gPDown = compProb(wl,L(i),qDown,Q(neighs),wneighs);
            
            grDown = (gPDown - pPrev(i))/stepSize;
            grDown = grDown * multMod;
                                    
            if (grDown > 0)
                
                if (q - grDown <= likelihoodLimit && q > likelihoodLimit)
                    if (dots(i) > .5)
                        c = 0;
                    else
                        QConnectivity = QNEW;
                        QConnectivity(i) = q - grDown;
                        c = isSimplePointMesh(fv,i,QConnectivity);
                    end
                else
                    c =1;
                end
                if (c == 1)
                    QNEW(i) = max(q - grDown, 0);
                    pPrev(i) = compProb(wl,L(i),QNEW(i),Q(neighs),wneighs);
                    downCounter = downCounter + 1;
                end
            else
                
                if (q - grDown > likelihoodLimit && q <= likelihoodLimit)
                    QConnectivity = QNEW;
                    QConnectivity(i) = q - grDown;
                    QConnectivity = 1+(QConnectivity *-1);
                    
                    c = isSimplePointMesh(fv,i,QConnectivity);
                else
                    c = 1;
                end
                if (c == 1)
                    QNEW(i) = min(q - grDown, 1);
                    pPrev(i) = compProb(wl,L(i),QNEW(i),Q(neighs),wneighs);
                    upCounter = upCounter + 1;
                end
                
            end
            
                        
        end
    end
    
    diff = sqrt(sum((Q-QNEW).^2));
    
    totalPr = sum(pPrev(L>0));
    totalPrDiffPerPoint = (totalPr - totalPr2)/sum(L>0);
    
    currFundPoints = sum(Q>0.5);
    
    if (totalPrDiffPerPoint < 0.0001 && (currFundPoints - prevFundPoints) == 0)
        endFlag = endFlag - 1;
    end
    totalPr2 = totalPr;
    disp([diff totalPrDiffPerPoint]);
    prevFundPoints = currFundPoints;
    
    Q = QNEW;
    disp(['Iteration: ' num2str(iterationCount) ' ** Up and down movement: ' num2str(upCounter) ', ' num2str(downCounter) ' ** Current fundus points (q > .5): ' num2str(currFundPoints)]);
end
