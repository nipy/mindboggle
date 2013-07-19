'''VALID'''

from numpy import zeros, array
from numpy.linalg import det

from .mon_comb_m import mon_comb

def D_CV_orig(facet,num_vertices,N,X,K,tri_matrix) :               #function [C,Vol] = D_CV_orig(facet,num_vertices,N,X,K,tri_matrix)

    C = zeros( (num_vertices,N+1,N+1,N+1) )                        #C=zeros(num_vertices,N+1,N+1,N+1);
                                                                   #
                                                                   #    % Vertex 1
                                                                   #    % Monomial combinations for each vertex in a given facet
    vertex_1 = X[K[facet,0],:] #Warn: 0-1                          #    vertex_1=X(K(facet,1),:);
                                                                   #    % Vertex
    vx = vertex_1[0] #Warn: 0-1                                    #    vx=vertex_1(1);
    vy = vertex_1[1] #Warn: 0-1                                    #    vy=vertex_1(2);
    vz = vertex_1[2] #Warn: 0-1                                    #    vz=vertex_1(3);
    tmp = (vx,vy,vz)
    C[0,:,:,:] = mon_comb(tri_matrix,tmp,N) # Warn: 0-1            #C(1,:,:,:)=mon_comb(tri_matrix,[vx vy vz],N);
                                                                   #
                                                                   #    % Vertex 2
                                                                   #    % Monomial combinations for each vertex in a given facet
    vertex_2 = X[K[facet,1],:] # Warn: 0-1                         #    vertex_2=X(K(facet,2),:);
                                                                   #    % Vertex
    vx = vertex_2[0] #Warn: 0-1                                    #    vx=vertex_2(1);
    vy = vertex_2[1] #Warn: 0-1                                    #    vy=vertex_2(2);
    vz = vertex_2[2] #Warn: 0-1                                    #    vz=vertex_2(3);
    tmp = (vx,vy,vz)
    C[1,:,:,:] = mon_comb(tri_matrix,tmp,N)                        #    C(2,:,:,:)=mon_comb(tri_matrix,[vx vy vz],N);
                                                                   #
                                                                   #    % Vertex 1
                                                                   #    % Monomial combinations for each vertex in a given facet
    vertex_3 = X[K[facet,2],:] #Warn: 0-1                          #    vertex_3=X(K(facet,3),:);
                                                                   #    % Vertex
    vx = vertex_3[0] #Warn: 0-1                                    #    vx=vertex_3(1);
    vy = vertex_3[1] #Warn: 0-1                                    #    vy=vertex_3(2);
    vz = vertex_3[2] #Warn: 0-1                                    #    vz=vertex_3(3);
    tmp = (vx,vy,vz)
    C[2,:,:,:] = mon_comb(tri_matrix,tmp,N) #Warn: 0-1             #    C(3,:,:,:)=mon_comb(tri_matrix,[vx vy vz],N);
                                                                   #    %Calculate Volume
    arr = array( [vertex_1, vertex_2, vertex_3] ).T
    print arr
    Vol = det(arr)    #    Vol=det(cat(2,vertex_1',vertex_2',vertex_3'));
    return C,Vol
