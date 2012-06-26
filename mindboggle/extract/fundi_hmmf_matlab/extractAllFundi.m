function fundi = extractAllFundi(fv,Depth,Cmean,Umin)

%inputs: 
% 1. fv: mesh represented by a struct having two fields: fv.faces and
% fv.vertices. fv.vertices has m x 3 elements, fv.faces has n x 3 elements.
% 
% 2. Depth: depth values in vector of size m x 1
% 
% 3. Cmean: mean curvature values in vector of size m x 1
%
% 4. Umin: directions of minimum curvature in matrix of size 3 x m
%
% output:
% 
% 1. fundi: a matrix of size m x indMax , where indMax is the size of
% sulci. each column represents a fundus for a specific sulcus. The values
% range between 0 and 1, with values above 0.5 considered as part of the
% fundus.

%depth threshold for defining sulci
depthTh = 0.2;

%extract sulci
Sulci = extractSulci(fv.faces,Depth,depthTh);
Sulci = fillSulciHoles(faces,Sulci);

%indMax is the number of sulci
indMax = max(Sulci(:));

% compute fundus likelihood values
disp('computing likelihood values');
L = zeros(length(Depth),1);
for ind = 1:indMax
    L = L + compFundusLikelihoodSingleSulcus(Cmean,Depth,Sulci,ind);
end

% extract fundi
fundi = zeros(length(L),indMax);
for ind = 1:indMax
    disp(ind);
    fundi(:,ind) = extractSingleFundus(L,Sulci,ind,fv,Umin);   
    save('fundiTemp20120624a.mat', 'fundi');   
end
