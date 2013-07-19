'''VALID'''

from numpy import zeros, power

def mon_comb(tri_matrix,vertex,N) :           #function c=mon_comb(tri_matrix,vertex,N)
                                              #% Computes the value for a specified monomial
                                              #% combination for a given vertex
                                              #
                                              #% Vertex coordinates
    x = vertex[0] #Warn: 0-1                  #x=vertex(1);
    y = vertex[1] #Warn: 0-1                  #y=vertex(2);
    z = vertex[2] #Warn: 0-1                                                     #z=vertex(3);
                                                                                 #% Computation
    N += 1 # Warn: hackery
    c = zeros((N,N,N)) #Warn: 0-1                                                #c=zeros(N+1,N+1,N+1);
    for i in xrange(N) :                                                         #for i=0:N
        for j in xrange(N-i) :                                                   #    for j=0:(N-i)
            for k in xrange(N-i-j) :                                             #        for k=0:(N-i-j)
                mon = power(x,i)*power(y,j)*power(z,k) #Need? 0-1    #            mon=power(x,i)*power(y,j)*power(z,k);
                c[i,j,k] = tri_matrix[i,j,k]*mon #Warn: 0-1                      #            c(i+1,j+1,k+1)=tri_matrix(i+1,j+1,k+1)*mon;
                                                                                 #        end
                                                                                 #    end
                                                                                 #end
    return c
