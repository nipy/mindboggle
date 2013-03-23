function a = learnDistributionParameters(feature,Sulci,LabelBoundarySulci,x)

NonLabelBoundarySulci = Sulci;
NonLabelBoundarySulci(LabelBoundarySulci == 1) = 0;

dataVec = feature(NonLabelBoundarySulci == 1);

if (min(x) < -0.9)
    inds0 = find(dataVec == 0);
    dataVec(inds0) = 0.05*randn(size(inds0,1),1);
end

[clMeansNL clSigmasNL WNL] = fitKNormalsOnHistogramEM(dataVec,x);


dataVec = feature(LabelBoundarySulci == 1);
if (min(x) < -0.9)
    inds0 = find(dataVec == 0);
    dataVec(inds0) = 0.05*randn(size(inds0,1),1);
end
[clMeansLB clSigmasLB WLB] = fitKNormalsOnHistogramEM(dataVec,x);

a = [clMeansNL clSigmasNL WNL clMeansLB clSigmasLB WLB];
