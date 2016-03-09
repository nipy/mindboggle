def D_SG_orig_part(cfg,num_facets,i,j,k,C,D,S,Vol,F) :                                  #function S = D_SG_orig_part(num_facets,i,j,k,C,D,S,Vol,F)
                                                                                        #
    for facet in cfg.rng(0,num_facets-1) : #!                                             #for facet=1:num_facets
                                                                                        #    % Sijk for a given facet
                                                                                        #    
                                                                                        #    
                                                                                        #    %aux_1 = factorial(i)*factorial(j)*factorial(k);
                                                                                        #    %aux_2 = factorial(i+j+k+2);
                                                                                        #    
        aux_1 = F[i]*F[j]*F[k] #!                                                       #    aux_1 = F(i + 1)*F(j + 1)*F(k + 1);
        aux_2 = F[i+j+k+2] #!                                                           #    aux_2 = F(i+j+k+2 + 1);
        aux = aux_1/aux_2                                                               #    aux=aux_1/aux_2;
                                                                                        #    
        tmp=0                                                                           #    tmp=0;
        for i1 in cfg.rng(0,i) :                                                        #    for i1=0:i
            for j1 in cfg.rng(0,j) :                                                    #        for j1=0:j
                for k1 in cfg.rng(0,k) :                                                #            for k1=0:k
                    tmp = tmp + (C[facet,i1,j1,k1]*D[facet,i-i1,j-j1,k-k1]) #!          #                tmp=tmp+(C(facet,i1+1,j1+1,k1+1)*D(facet,i-i1+1,j-j1+1,k-k1+1));
                                                                                        #            end
                                                                                        #        end
                                                                                        #    end
        tmp = tmp*aux                                                                   #    tmp=tmp*aux;
                                                                                        #    
        S[facet,i,j,k]=(Vol[facet]*tmp)/(i+j+k+3) #!                                    #    S(facet,i+1,j+1,k+1)=(Vol(facet)*tmp)/(i+j+k+3);
                                                                                        #end
    return S
