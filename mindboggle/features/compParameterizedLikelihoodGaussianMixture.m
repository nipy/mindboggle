function L = compParameterizedLikelihoodGaussianMixture(aDepth,aCurv,Sulci,norm_depths,curvature)

L = zeros(length(Sulci),1);
ProbNL = L;
ProbLB = L;

k = 2;
normNL = 1/(((2*pi)^(k/2)) .* aDepth(:,2).*aCurv(:,2));
normLB = 1/(((2*pi)^(k/2)) .* aDepth(:,5).*aCurv(:,5));

nmix = size(aDepth,1);

for j = 1:nmix
    for i = 1:length(Sulci)
        if(Sulci(i) > -1)
            NLexp = aDepth(j,3) * ((norm_depths(i)-aDepth(j,1))^2)/(aDepth(j,2)^2);
            NLexp = NLexp + aCurv(j,3) * ((curvature(i)-aCurv(j,1))^2)/(aCurv(j,2)^2);
            NLexp = NLexp * (-1/2);
            ProbNL(i) = ProbNL(i) + normNL(j) * exp(NLexp);
            
            LBexp = aDepth(j,6) * ((norm_depths(i)-aDepth(j,4))^2)/(aDepth(j,5)^2);
            LBexp = LBexp + aCurv(j,6) * ((curvature(i)-aCurv(j,4))^2)/(aCurv(j,5)^2);
            LBexp = LBexp * (-1/2);
            ProbLB(i) = ProbLB(i) + normLB(j) * exp(LBexp);
        end
    end
end

L = ProbLB ./ (ProbNL + ProbLB);

L(norm_depths < 0) = 0;
