function Q = extractPitsHMMF(fv,fvInflated,Depth,pitInds,datapath,lambda,lh)

geodesic = 0;

if (length(datapath) < 2)    
    saveResults = 0;
else
    saveResults = 1;
end
tic


if (pitInds(1) < 1)
    
    %find local maxima of depth as pit candidates
    disp('Finding local maxima of depth');
    %[vertices faces curv pitInds pits] = xtractPitsSurf(fv.faces,fv.vertices,Depth);
    [vertices faces pitInds] = findPitCandidates(fv.faces,fv.vertices,Depth);
end

%compute curvature maps
disp('Computing curvature');
options.curvature_smoothing = 8;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(fv.vertices,fv.faces,options);

% compute pit likelihood for each point
disp('Computing pit likelihood');
L = compPitLikelihood(Cgauss,Cmean,Depth);

pits = fvInflated.vertices(pitInds,:);


% %%%%%%%%%%
% disp('Likelihoods:');
% disp(L(pitInds));
% disp(pits);
%%%%%%%%%%

if (geodesic)
    
    % compute weight matrix for geodesic distance computation
    disp('Computing weight matrix');
    W = compWeightMatrix(fv);
    
    % compute geodesic distance to nearby candidates
    disp('Computing geodesic distances');
    D = compGeodesicDistanceToMaxDist(pitInds,pits,W);
    
else
    D = compEuclidianDistanceForPointGroup(pitInds,pits);
end

% compute optimal measure field
disp('Computing optimal label field');
Q = compOptimalHMMF(L,D,pitInds,lambda);

if (saveResults)
    if (lh == 1) 
        save([datapath '/miccaiData/pitResultsLH' num2str(1000*lambda)]);
    else
        save([datapath '/miccaiData/pitResultsRH' num2str(1000*lambda)]);
    end
end

toc
