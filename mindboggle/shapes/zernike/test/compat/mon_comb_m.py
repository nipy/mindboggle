def mon_comb(cfg,tri_matrix,vertex,N) :                          #function c=mon_comb(tri_matrix,vertex,N)
                                                                 #% Computes the value for a specified monomial
                                                                 #% combination for a given vertex
                                                                 #
                                                                 #% Vertex coordinates
    x = vertex[0] #!                                             #x=vertex(1);
    y = vertex[1] #!                                             #y=vertex(2);
    z = vertex[2] #!                                             #z=vertex(3);
                                                                 #% Computation
    c=cfg.zeros(N+1,N+1,N+1)                                     #c=zeros(N+1,N+1,N+1);
    for i in cfg.rng(0,N) :                                      #for i=0:N
        for j in cfg.rng(0,N-i) :                                #    for j=0:(N-i)
            for k in cfg.rng(0,N-i-j) :                          #        for k=0:(N-i-j)
                mon=cfg.power(x,i)*cfg.power(y,j)*cfg.power(z,k) #            mon=power(x,i)*power(y,j)*power(z,k);
                c[i,j,k] = tri_matrix[i,j,k]*mon #!              #            c(i+1,j+1,k+1)=tri_matrix(i+1,j+1,k+1)*mon;
                                                                 #        end
                                                                 #    end
                                                                 #end
    return c
