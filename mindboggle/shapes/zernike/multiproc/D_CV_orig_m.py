def D_CV_orig(cfg,facet,num_vertices,N,X,K,tri_matrix) :      #function [C,Vol] = D_CV_orig(facet,num_vertices,N,X,K,tri_matrix)
                                                              #
    C = cfg.zeros(num_vertices,N+1,N+1,N+1)                   #C=zeros(num_vertices,N+1,N+1,N+1);
                                                              #
                                                              #    % Vertex 1
                                                              #    % Monomial combinations for each vertex in a given facet
    vertex_1=X[K[facet-1,0]-1,:] #! facet is 1-indexed        #    vertex_1=X(K(facet,1),:);
                                                              #    % Vertex
    vx = vertex_1[0] #!                                       #    vx=vertex_1(1);
    vy = vertex_1[1] #!                                       #    vy=vertex_1(2);
    vz = vertex_1[2] #!                                       #    vz=vertex_1(3);
    C[0,:,:,:]=cfg.mon_comb(tri_matrix,[vx,vy,vz],N) #!       #    C(1,:,:,:)=mon_comb(tri_matrix,[vx vy vz],N);
                                                              #
                                                              #    % Vertex 2
                                                              #    % Monomial combinations for each vertex in a given facet
    vertex_2=X[K[facet-1,1]-1,:] #!                           #    vertex_2=X(K(facet,2),:);
                                                              #    % Vertex
    vx=vertex_2[0] #!                                         #    vx=vertex_2(1);
    vy=vertex_2[1] #!                                         #    vy=vertex_2(2);
    vz=vertex_2[2] #!                                         #    vz=vertex_2(3);
    C[1,:,:,:]=cfg.mon_comb(tri_matrix,[vx,vy,vz],N) #!       #    C(2,:,:,:)=mon_comb(tri_matrix,[vx vy vz],N);
                                                              #
                                                              #    % Vertex 1
                                                              #    % Monomial combinations for each vertex in a given facet
    vertex_3=X[K[facet-1,2]-1,:] #!                           #    vertex_3=X(K(facet,3),:);
                                                              #    % Vertex
    vx=vertex_3[0] #!                                         #    vx=vertex_3(1);
    vy=vertex_3[1] #!                                         #    vy=vertex_3(2);
    vz=vertex_3[2] #!                                         #    vz=vertex_3(3);
    C[2,:,:,:]=cfg.mon_comb(tri_matrix,[vx,vy,vz],N) #!       #    C(3,:,:,:)=mon_comb(tri_matrix,[vx vy vz],N);
                                                              #    %Calculate Volume
    Vol=cfg.det(cfg.cat(1,vertex_1,vertex_2,vertex_3)) #!axis #    Vol=det(cat(2,vertex_1',vertex_2',vertex_3'));
    return C,Vol
