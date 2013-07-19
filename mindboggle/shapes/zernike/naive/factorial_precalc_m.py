'''VALID'''

def factorial_precalc(N) :                                 #function F = factorial_precalc(N)
                                                           #
    if N == 0 : return zeros(N,dtype='int')
    F = zeros(N+3,dtype='float')
    for i in xrange(N) :                                   #for i = 1:N    
        for j in xrange(N-i) :                             #    for j=0:(N-i)
            for k in xrange(N-i-j) :                       #        for k=0:(N-i-j)
                                                           #
                ii = i+1
                F[ii] = factorial(ii)                      #            F(i + 1) = factorial(i);
                F[j] = factorial(j)                        #            F(j + 1) = factorial(j);
                F[k] = factorial(k)                        #            F(k + 1) = factorial(k);
                F[ii+j+k+2] = factorial(ii+j+k+2)          #            F(i+j+k+2 + 1) = factorial(i+j+k+2);
                                                           #            
                                                           #        end
                                                           #    end
                                                           #end
    return F
