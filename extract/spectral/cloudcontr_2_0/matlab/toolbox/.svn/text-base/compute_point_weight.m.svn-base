function W = compute_point_weight(pts, type, rings, options)

% compute_point_weight - compute a weight matrix
%
%   W = compute_point_weight(pts, type, rings, options);
%
%   W is sparse weight matrix and W(i,j)=0 is vertex i and vertex j are not
%   connected.
%
%   type is either 
%       'combinatorial': W(i,j)=1 is vertex i is conntected to vertex j.
%       'distance': W(i,j) = 1/d_ij^2 where d_ij is distance between vertex
%           i and j.
%       'spring': W(i,j) = 1/d_ij where d_ij is distance between vertex
%           i and j.
%       'conformal' or 'dcp': W(i,j) = cot(alpha_ij)+cot(beta_ij) where alpha_ij and
%           beta_ij are the adjacent angle to edge (i,j). Refer to Skeleton
%           Extraction by Mesh Extraction_08, and Intrinsic Parameterizations of Surface Meshes_02.
%       'Laplace-Beltrami': W(i,j) = (cot(alpha_ij)+cot(beta_ij))/2 where alpha_ij and
%           beta_ij are the adjacent angle to edge (i,j). Refer to Computing discrete minimal surfaces
%           andtheir conjugates_93, Lemma 2 of On the convergence of metric and geometric properties of
%           polyhedral surfaces_06, and Characterizing Shape Using
%           Conformal Factors_08, and ,
%       'mvc': W(i,j) = [tan(/_kij/2)+tan(/_jil/2)]/d_ij where /_kij and /_jil
%           are angles at i
%
%   If options.normalize=1, the the rows of W are normalize to sum to 1.
%   If M.ring is offered, we can avoid compute it.
%   
%   NB: we adapt it for better numeric stablity for contracting a model by adding the following lines
%     tmp = sum(W(i,:));
%     if tmp>10000
%         W(i,:) = W(i,:)*10000/tmp;
%     end
%   Copyright (c) 2009 jjcao
options.null = 0;

if nargin<2
    type = 'conformal';
end
        
switch lower(type)
    case 'combinatorial'
        W = compute_point_weight_combinatorial(pts, rings);
    case 'distance'
       warning('not implemented!'); 
    case 'spring'
        W = compute_point_weight_spring(pts, rings);        
    case {'conformal','dcp'} % conformal laplacian  
        W = compute_point_weight_dcp(pts, rings);
    case 'laplace-beltrami' %
        W = compute_point_weight_dcp(pts, rings)*0.5;
    case 'mvc'% mvc laplacian
        W = compute_point_weight_mvc(pts, rings);  
    otherwise
        error('Unknown type!!')
end

%#########################################################################
function W = compute_point_weight_combinatorial(points, rings)
n = length(points);
W = sparse(n,n);
for i = 1:n
    ring = rings{i};
    if ring(1) == ring(end)
        ring = ring(1,1:(end-1));
    end
    for j = ring
        W(i,j) = 1.0;
    end
end
function W = compute_point_weight_spring(points, rings)
n = length(points);
W = sparse(n,n);
for i = 1:n
    vi = points(i,:);
    ring = rings{i};
    if ring(1) == ring(end)
        ring = ring(1,1:(end-1));
    end
    for j = ring                
        vj = points(j,:);        
        W(i,j) = 1./sqrt(sum((vi-vj).^2));
    end
    tmp = sum(W(i,:));
    if tmp>10000
        W(i,:) = W(i,:)*10000/tmp;
    end
end
%#########################################################################
function W = compute_point_weight_dcp(points, rings)
n = length(points);
W = sparse(n,n);
for i = 1:n
    ring = rings{i};

    tmp = size(ring,2)-1;
    for ii = 1: tmp
        j = ring(ii); k = ring(ii+1);
        vi = points(i,:);
        vj = points(j,:);
        vk = points(k,:);
        
        % new % Oscar08 use this 
        u = vk-vi; v = vk-vj;
        cot1 = dot(u,v)/norm(cross(u,v));
        W(i,j) = W(i,j) + cot1;
        u = vj-vi; v = vj-vk;
        cot2 = dot(u,v)/norm(cross(u,v));
        W(i,k) = W(i,k) + cot2;        
%         % old
%         % angles
%         alpha = myangle(vk-vi,vk-vj);
%         beta = myangle(vj-vi,vj-vk);
%         % add weight
%         W(i,j) = W(i,j) + cot( alpha );
%         W(i,k) = W(i,k) + cot( beta );
    end
    
    tmp = abs(sum(W(i,:)));
    if tmp>10000
        W(i,:) = W(i,:)*10000/tmp;
    end
end
function W = compute_point_weight_mvc(points, rings)
n = length(points);
W = sparse(n,n);
for i = 1:n
    ring = rings{i};

    tmp = size(ring,2)-1;
    for ii = 1: tmp
        j = ring(ii); k = ring(ii+1);
        vi = points(i,:);
        vj = points(j,:);
        vk = points(k,:);
        
        % angles
        alpha = myangle(vi-vk,vi-vj);
        % add weight
        W(i,j) = W(i,j) + tan( 0.5*alpha )/sqrt(sum((vi-vj).^2));
        W(i,k) = W(i,k) + tan( 0.5*alpha )/sqrt(sum((vi-vk).^2));
    end
    
    tmp = sum(W(i,:));
    if tmp>10000
        W(i,:) = W(i,:)*10000/tmp;
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = myangle(u,v);

du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) ); 