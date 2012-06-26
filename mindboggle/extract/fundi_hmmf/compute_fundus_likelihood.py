function L = compFundusLikelihoodSingleSulcus(Cmean,Depth,Sulci,ind)

massThr1 = .6;
massThr2 = .05;
massThr3 = .3;
highMapValue = .9;

Sulci(Sulci ~= ind) = 0;

%switch curvature values to opposite
Cmean = (Cmean*-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

depthSearch = 0;
massLeft = 1;

% find depth value where less than [massThr1] are larger
while(massLeft > massThr1)
    depthSearch = depthSearch + .01;
    massLeft = sum(Depth(Sulci > 0) > depthSearch) / sum(Sulci > 0);
end
halfValDepth = depthSearch;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find depth value where less than [massThr2] are larger
while(massLeft > massThr2)
    depthSearch = depthSearch + .01;
    massLeft = sum(Depth(Sulci > 0) > depthSearch) / sum(Sulci > 0);
end
slopeDepth = (-1/(depthSearch - halfValDepth)) * log((1/highMapValue)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find slope similarly for curvature values
massLeft = 1;
depthSearch = 0;
while(massLeft > massThr3)
    depthSearch = depthSearch + .0001;    
    massLeft = sum(Cmean(Sulci > 0) > depthSearch) / sum(Sulci > 0);
end
slopeCmean = (-1/(depthSearch)) * log((1/highMapValue)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Map curvature and depth values with sigmoidal functions to range [0,1]
st2 = 1./(1+exp(-slopeCmean*Cmean));
st3 = 1./(1+exp(-slopeDepth*(Depth - halfValDepth)));

L = st2 .* st3;
L(Sulci ==0) = 0;
