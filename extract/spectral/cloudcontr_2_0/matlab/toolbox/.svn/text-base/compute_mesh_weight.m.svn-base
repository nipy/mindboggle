function W = compute_mesh_weight(vertices,faces,type,options)

% compute_mesh_weight - compute a weight matrix
%
%   W = compute_mesh_weight(vertices,faces,type,options);
%
%   W is sparse weight matrix and W(i,j)=0 is vertex i and vertex j are not
%   connected in the mesh.
%
%   todo:
%       1. validate spring weights
%       2. add Mean_curvature weights
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
%           polyhedral surfaces_06, and Characterizing Shape Using Conformal Factors_08, and ,
%       'Mean_curvature',W(i,j) = (1/area_i)*(cot(alpha_ij)+cot(beta_ij))/2 where alpha_ij and
%           beta_ij are the adjacent angle to edge (i,j), area_i is the
%           area of vertex i's Voroni vicinity. Refer to ??
%       'mvc': W(i,j) = [tan(/_kij/2)+tan(/_jil/2)]/d_ij where /_kij and /_jil
%           are angles at i
%
%   If options.ring is offered, the computation of it can be avoided.
%
%   Add spring and mvc weight (JJCAO, 2009)
%   Add options.ring (JJCAO, 2009)
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
[vertices,faces] = check_face_vertex(vertices,faces);

n = max(max(faces));

if isfield(options, 'verb')
    verb = options.verb;
else
    verb = n>5000;
end

if nargin<3
    type = 'dcp';
end

switch lower(type)
    case 'combinatorial'
        W = triangulation2adjacency(faces);
    case 'distance'
        W = my_euclidean_distance_2(triangulation2adjacency(faces),vertices);
        W(W>0) = 1./W(W>0);
    case 'spring'
        W = sqrt(my_euclidean_distance(triangulation2adjacency(faces),vertices));
        W(W>0) = 1./W(W>0);
    case {'conformal','dcp'} % conformal laplacian  
        W = compute_mesh_weight_dcp(vertices, faces, verb);        
    case 'laplace-beltrami' %
        W = compute_mesh_weight_dcp(vertices, faces, verb)*0.5;            
    case 'mvc'% mvc laplacian
        if isfield(options, 'rings')
            rings = options.rings;
        else
            rings = compute_vertex_face_ring(faces);
        end
        W = sparse(n,n);
        for i = 1:n
            if verb
                progressbar(i,n);
            end
            for b = rings{i}
                % b is a face adjacent to a
                bf = faces(:,b);
                % compute complementary vertices
                if bf(1)==i
                    v = bf(2:3);
                elseif bf(2)==i
                    v = bf([1 3]);
                elseif bf(3)==i
                    v = bf(1:2);
                else
                    error('Problem in face ring.');
                end
                j = v(1); k = v(2);
                vi = vertices(:,i);
                vj = vertices(:,j);
                vk = vertices(:,k);
                % angles
                alpha = myangle(vi-vk,vi-vj);
                % add weight
                W(i,j) = W(i,j) + tan( 0.5*alpha )/sqrt(sum((vi-vj).^2));
                W(i,k) = W(i,k) + tan( 0.5*alpha )/sqrt(sum((vi-vk).^2));
            end
        end     
    otherwise
        error('Unknown type.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = compute_mesh_weight_dcp(verts, faces, verb)
n = size(verts,2);
W = sparse(n,n);
% new implementation
for i = 1:size(faces,2)
    c = faces(:,i);
    
    v1 = verts(:,c(1));
    v2 = verts(:,c(2));
    v3 = verts(:,c(3));
    
    cot1 = dot(v2-v1,v3-v1)/norm(cross(v2-v1,v3-v1));
    cot2 = dot(v3-v2,v1-v2)/norm(cross(v3-v2,v1-v2));
    cot3 = dot(v1-v3,v2-v3)/norm(cross(v1-v3,v2-v3));
    
    W(c(2),c(3)) = W(c(2),c(3)) + cot1;
    W(c(3),c(2)) = W(c(3),c(2)) + cot1;
    W(c(3),c(1)) = W(c(3),c(1)) + cot2;
    W(c(1),c(3)) = W(c(1),c(3)) + cot2;
    W(c(1),c(2)) = W(c(1),c(2)) + cot3;
    W(c(2),c(1)) = W(c(2),c(1)) + cot3;    
end

%%
% old implementation
% for i = 1:n
%     if verb
%         progressbar(i,n);
%     end
%     for b = rings{i}
%         % b is a face adjacent to a
%         bf = faces(:,b);
%         % compute complementary vertices
%         if bf(1)==i
%             v = bf(2:3);
%         elseif bf(2)==i
%             v = bf([1 3]);
%         elseif bf(3)==i
%             v = bf(1:2);
%         else
%             error('Problem in face ring.');
%         end
%         j = v(1); k = v(2);
%         vi = verts(:,i);
%         vj = verts(:,j);
%         vk = verts(:,k);
%         
%         u = vk-vi; v = vk-vj;
%         W(i,j) = W(i,j) + dot(u,v)/norm(cross(u,v));
%         u = vj-vi; v = vj-vk;
%         W(i,k) = W(i,k) + dot(u,v)/norm(cross(u,v));
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = myangle(u,v)
du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = my_euclidean_distance_2(A,vertex)
% square euclidean distance
if size(vertex,1)<size(vertex,2)
    vertex = vertex';
end

[i,j,s] = find(sparse(A));
d = sum( (vertex(i,:) - vertex(j,:)).^2, 2);
W = sparse(i,j,d);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = my_euclidean_distance(A,vertex)
% euclidean distance
if size(vertex,1)<size(vertex,2)
    vertex = vertex';
end

[i,j,s] = find(sparse(A));
d = sum( (vertex(i,:) - vertex(j,:)).^2, 2).^0.5;
W = sparse(i,j,d);  