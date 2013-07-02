'''VALID'''

def trinomial_matrix(N) :                                            #function tri_matrix=trinomial_matrix(N)
                                                                     #% Computes the trinomial of
                                                                     #% the three input arguments
    
    N += 1 #Warn hackery
    tri_matrix = numpy.zeros((N,N,N),dtype='int') #Warn: 0-1         #tri_matrix=zeros(N+1,N+1,N+1);
    for i in xrange(N) :                                             #for i=0:N
        for j in xrange(N-i) :                                       #for j=0:(N-i)
            for k in xrange(N-i-j) :                                 #for k=0:(N-i-j)
                tri_matrix[i,j,k] = trinomial(i,j,k) #Warn: 0-1      #tri_matrix(i+1,j+1,k+1)=trinomial(i,j,k);
                                                                     #end
                                                                     #end
                                                                     #end
    return tri_matrix
