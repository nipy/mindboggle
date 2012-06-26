function fundus = extractSingleFundus(L,Sulci,ind,fv,Umin)

% determine minimum size for a sulcus from which we want to find fundi.
% This can be set to 0, the result will just have very short fundi.
minSulcusSize = 30;

L2 = zeros(size(L));
L2(Sulci==ind) = L(Sulci==ind);
clear L;

sizeCheck = sum(L2>0.5);

if (sizeCheck > minSulcusSize)
    tic;
    Q = findSupportPoints(fv,L2,Umin);
    %%%%%%%
    nList = compNeighborList(fv.faces,find(L2>0));
    shortInds = find(L2>0);
    %%%%%%%%%%%%%%%
    %initialize all likelihoodvalues within sulcus between [.5 1.0]. This
    %is necessary to guarantee correct topology
    initValues = (L2 + 1.001)/2;
    initValues(L2 == 0) = 0;
    initValues(initValues > 1) = 1;
    %%%%%%%
    disp(find(Q>0.5));
    %Q = connectTheDotsV2(L2,initValues,fv,Q,nList,shortInds);
    if (sum(Q>0.5) > 1)
        Q = connectTheDotsV3(L2,initValues,fv,Q,nList,shortInds);        
        fundus = Q;
    else
        fundus = zeros(size(Q));
    end
    toc;
else
    disp('Not Enough Values');
    disp(sizeCheck);
    fundus = zeros(size(L2));
end
