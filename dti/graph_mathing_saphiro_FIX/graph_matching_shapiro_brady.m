function [associations,Z] = graph_matching_shapiro_brady(adj_1,adj_2,options)
% This function implement the algorithm for the wieghted graph matching
% proposed in:
%
% Larry S Shapiro and J Michael Brady "Feature-based correspondence: an
% eigenvector approach" Image and vision computing, 10 (1992) 283-288.
%
% NB: elements in the adjacency matrices are interpreted as distances
% between nodes.
% Ex: adj_1(r,c) is the distance between nodes r and c in the graph 1.

% Options:
% options.alignment_method : to select the metod to align the eigenvectors
%                            between the two bases.
%                            0: no alignment is performed
%                            1: the absolute value for each element is
%                               computed
%                            2: the correlation of the progressive area on
%                               the element frequency is used to decide the
%                               alignment versus.

%                                                       1) Check input data
% -------------------------------------------------------------------------
[nR1,nC1] = size(adj_1);
if nR1~=nC1
    disp('ERROR: adjacency matrices are supposed to be squared')
end
[nR2,nC2] = size(adj_2);
if nR2~=nC2
    disp('ERROR: adjacency matrices are supposed to be squared')
end
final_nC = min([nC1,nC2]);

%                                       2) Computing the proximity matrices
% -------------------------------------------------------------------------
H1 = adj2prox(adj_1,options);
H2 = adj2prox(adj_2,options);

%                                3) Decomposition of the proximity matrices
% -------------------------------------------------------------------------
[V1,D1] = eig_decompose(H1);
[V2,D2] = eig_decompose(H2);

% Truncating the least significant modes in order to achieve the same
% amount of columns between V_1 and V_2
V1_t = V1(:,1:final_nC);
V2_t = V2(:,1:final_nC);

% Eigenvalues base must be aligned
[V1_ta,V2_ta] = alignEigenvetors(V1_t,V2_t,options);

%                                       4) Computing the association matrix
% -------------------------------------------------------------------------
Z = get_association_matriz(V1_ta,V2_ta);
associations = compute_associations(Z);





% -------------------------------------------------------------------------
% |                            Other Functions                            |
% -------------------------------------------------------------------------
function [H] = adj2prox(adj,options)
% This function compute the proximity matrix given the adjacency matrix
% using the Gaussian-weighted metric proposed in the paper.

[nR,nC]=size(adj);
%sigma = options.sigma;
sigma = mean(adj(find(adj)));
H = zeros(nR,nC);
for r=1:nR
    for c=1:nC
        if r~=c 
            H(r,c) = exp(-(adj(r,c)^2)/(2*sigma));
        else
            % Element in the main diagonal are all 1, which indicates the
            % maximum proximity.
            H(r,c)=1;
        end
    end
end





function [V,D] = eig_decompose(H)
% This function decompose the given proximity matrix using the
% eigenvectors. In particular column of V are the eigenvectors of H ordered
% using the eigenvalues. D is a diagonal matrix containing the eigenvalues,
% from the largest one to the smallest one.
% NOTE: H = V*D*V'

[V,D]=eig(H);
% Note: matlab gives the eignevalues from the smallest to the largest.
V = fliplr(V);
D = fliplr(flipud(D));





function [V1a,V2a] = alignEigenvetors(V1,V2,options)
% The two basis of eigenvectors must be aligned in order to compute the
% association matrix
% NOTE: this part is not well explained in the paper.
version_to_use = options.alignment_method;
switch version_to_use
    case 0
        % Version 0: doing nothing
        V1a = V1;
        V2a = V2;
    case 1
        % Version 1: I don't know how align the two basis, however later
        % only the puntual difference is used, therefore the absolute
        % values of V1 and V2 can also be computed. The drawback is an
        % increased sensitivity to noise.
        V1a = abs(V1);
        V2a = abs(V2);
    case 2
        % Version 2: the idea is to check sign of the second eigenvector
        % that better match the distribution of values of the first
        % eigenvector.
        % In this method each eigenvetor of V2 is aligned to the
        % coresponding eigenvector of V1 independently using the
        % correlation on the histogram area vector.
        nVett = size(V2,2);
        V1a = V1;
        V2a = zeros(size(V2));
        xVett = [min([-1, min(min(V2)), min(min(V1))]):0.05:max([-1, max(max(V2)), max(max(V1))])]';
        len_xVett = length(xVett);
        for col =1:nVett
            vett1  = V1(:,col);
            vett2p = V2(:,col); % positive version of the V2 column
            vett2n = V2(:,col); % negative version of the V2 column
            hist1  = hist(vett1,xVett)./length(vett1); % computing the pseudo frequency of the elements in vett1
            hist2p = hist(vett2p,xVett)./length(vett2p); % computing the pseudo frequency of the elements in vett2p
            hist2n = hist(vett2n,xVett)./length(vett2n); % computing the pseudo frequency of the elements in vett2p            
            areaVett1  = zeros(len_xVett,1);
            areaVett2p = zeros(len_xVett,1);
            areaVett2n = zeros(len_xVett,1);
            for cont = 2:len_xVett
                areaVett1(cont) = trapz(xVett(1:cont),hist1(1:cont));
                areaVett2p(cont) = trapz(xVett(1:cont),hist2p(1:cont));
                areaVett2n(cont) = trapz(xVett(1:cont),hist2n(1:cont)); 
            end
            corrMat = corrcoef(areaVett1,areaVett2p);
            positive_coeff = corrMat(2,1);
            corrMat = corrcoef(areaVett1,areaVett2n);
            negative_coeff = corrMat(2,1);
            if positive_coeff>negative_coeff
                V2a(:,col) = vett2p;
            else
                V2a(:,col) = vett2n;
            end
        end
end





function [Z] = get_association_matriz(V1,V2)
% Given the matrices V1 and V2 already truncated and aligned the method
% compute the association matrix Z following the method suggested in the
% paper.

[nR1,nC1]=size(V1);
[nR2,nC2]=size(V2);
if nC1~=nC2
    % Matrices V1 and V2 must have the same number of columns
    disp('ERROR: to compute the association matrix the matrices V1 and V2 must have the same number of columns')
    Z=0;
else
    Z = zeros(nR1,nR2);
    for r1=1:nR1
        for r2 = 1:nR2
            Z(r1,r2) = norm(V1(r1,:)-V2(r2,:)).^2;
        end
    end
end





function [associations] = compute_associations(Z)
% Following the paper, Z contains an association measure for each
% couple of nodes between the two graph. An association between nodes i,j
% is found is Z(i,j) is the minimum values along i-th row and j-th column.

associations = [];
pos =1;
for r =1:size(Z,1)
    [minValueR,minPosR]=min(Z(r,:));
    [minValueC,minPosC]=min(Z(:,minPosR));
    if minPosC==r
        % An association is found
        associations(pos,1) = r;
        associations(pos,2) = minPosR;
        pos = pos+1;
    end
end