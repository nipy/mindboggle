def trinomial_matrix(cfg,N) :                                  #function tri_matrix=trinomial_matrix(N)
                                                               #% Computes the trinomial of
                                                               #% the three input arguments
                                                               #
    tri_matrix=cfg.zeros(N+1,N+1,N+1)                          #tri_matrix=zeros(N+1,N+1,N+1);
    for i in cfg.rng(0,N) :                                    #for i=0:N
        for j in cfg.rng(0,N-i) :                              #    for j=0:(N-i)
            for k in cfg.rng(0,N-i-j) :                        #        for k=0:(N-i-j)
                tri_matrix[i,j,k] = cfg.trinomial(i,j,k) #!    #            tri_matrix(i+1,j+1,k+1)=trinomial(i,j,k);
                                                               #        end
                                                               #    end
                                                               #end
    return tri_matrix
