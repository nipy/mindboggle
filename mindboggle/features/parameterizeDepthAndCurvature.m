function [aDepth, aCurv] = parameterizeDepthAndCurvature(norm_depths,curvature,label_borders_sulci,Sulci)

x = 0:0.02:1;
aDepth = learnDistributionParameters(norm_depths,Sulci,label_borders_sulci,x);

x = -1:0.02:1;
aCurv = learnDistributionParameters(curvature,Sulci,label_borders_sulci,x);
