'''VALID'''

def D_SG_orig(num_facets,i,N,C,D,Vol,F) :                           #function G = D_SG_orig(num_facets,i,N,C,D,Vol,F)
                                                                    #
    G = zeros((N+1,N+1))                                            #G=zeros(N+1,N+1);
    S = zeros((num_facets,N+1,N+1,N+1))                             #S=zeros(num_facets,N+1,N+1,N+1);
                                                                    #
    for j in xrange(N-i+1) : #Warn: xrange                          #for j=0:(N-i)
        for k in xrange(N-i-j+1) : #Warn: xrange                    #    for k=0:(N-i-j)
            print j,k
                                                                    #        
            S = D_SG_orig_part(num_facets,i,j,k,C,D,S,Vol,F)        #        S = D_SG_orig_part(num_facets,i,j,k,C,D,S,Vol,F);
                                                                    #        % Geometric moments after summing through all facets
            G[j,k] = sum(S[:,i,j,k]) # Warn: 0-1                    #        G(j+1,k+1)=sum(S(:,i+1,j+1,k+1));
                                                                    #    end
                                                                    #end
    return G

