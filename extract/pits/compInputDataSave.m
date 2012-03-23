function a = compInputDataSave(fv,fvInflated,Depth,datapath,lh)

maxDist = 2.0;

if (length(datapath) < 2)    
    saveResults = 0;
else
    saveResults = 1;
end

%compute curvature maps
disp('Computing curvature');
options.curvature_smoothing = 8;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(fv.vertices,fv.faces,options);

% compute pit likelihood for each point
disp('Computing pit likelihood');
L = compPitLikelihood(Cgauss,Cmean,Depth);

[pitInds pits] = findPitCandidatesLocalMaxWithMaxDistLikelihoodLimit(fvInflated.faces,fvInflated.vertices,Depth,L,maxDist);
%pits = fvInflated.vertices(pitInds,:);

%D = compEuclidianDistanceForPointGroup(pitInds,pits);

if (saveResults)
    if (lh == 1) 
        save([datapath '/miccaiData/inputv6LH']);
    else
        save([datapath '/miccaiData/inputv6RH']);
    end
end

a = 0;
