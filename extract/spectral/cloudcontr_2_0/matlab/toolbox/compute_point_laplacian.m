function L = compute_point_laplacian(pts, type, rings, options)

% compute_point_laplacian - compute a laplacian matrix L, L = -(Laplacian operator).
%
%   L = compute_point_laplacian(pts, type, rings, options)
%
%   If options.symmetrize=1 and options.normalize=0 then 
%       L = D-W
%   If options.symmetrize=1 and options.normalize=1 then 
%       L = eye(n)-D^{-1/2}*W*D^{-1/2}
%   If options.symmetrize=0 and options.normalize=1 then 
%       L = eye(n)-D^{-1}*W.
%   where D=diag(sum(W,2)) and W is the unormalized weight matrix 
%   (see compute_mesh_weight).
%
%   type can 'combinatorial', 'distance', 'conformal'.
%
%   See also compute_point_weight.
%
%   Changed default value for symmetrize from 1 to 0 (JJCAO)
%   Copyright (c) 2009 jjcao

options.null = 0;
normalize = getoptions(options, 'normalize', 0);
symmetrize = getoptions(options, 'symmetrize', 1);

W = compute_point_weight(pts, type, rings, options);
n = size(W,1);
    
if symmetrize==1 && normalize==0
    L = diag(sum(W,2)) - W;
elseif symmetrize==1 && normalize==1
    L = speye(n) - diag(sum(W,2).^(-1/2)) * W * diag(sum(W,2).^(-1/2));
elseif symmetrize==0 && normalize==1
    L = speye(n) - diag(sum(W,2).^(-1)) * W;
else
    error('Does not work with symmetrize=0 and normalize=0');    
end
