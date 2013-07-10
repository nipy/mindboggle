import profilehooks
import numpy
import sys
import time
import itertools

#@profilehooks.timecall
def itertools_product(cfg,C,N) :
    D = cfg.zeros(N+1,N+1,N+1)                                                          #D = zeros(N+1,N+1,N+1);
    for a,b,c in cfg.irng((0,N),repeat=3 ) :                                            #for a=0:N
                                                                                        #    for b=0:N
                                                                                        #        for c=0:N
        temp = 0                                                                        #            temp=0;
        for i2,j2,k2 in cfg.irng((0,a),(0,b),(0,c)) :                                   #            for i2=0:a
                                                                                        #                for j2=0:b
                                                                                        #                    for k2=0:c
            temp += C[1,i2,j2,k2]*C[2,a-i2,b-j2,c-k2] #!                                #                        temp=temp+C(2,i2+1,j2+1,k2+1)*C(3,a-i2+1,b-j2+1,c-k2+1);
                                                                                        #                    end
                                                                                        #                end
                                                                                        #            end          
        D[a,b,c] = temp #!                                                              #            D(a+1,b+1,c+1)=temp;
                                                                                        #        end
                                                                                        #    end
                                                                                        #end
    return D


#@profilehooks.timecall
def orig(cfg,C,N) :                                                                #function D=Dabc_orig(C,N)
                                                                                        #% Pre-computes the values for D
                                                                                        #% and store them in matrix form
                                                                                        #% to be used for Geometric Moments
                                                                                        #
    #f = file('orig.indx','a')
    D = cfg.zeros(N+1,N+1,N+1)                                                          #D = zeros(N+1,N+1,N+1);
    for a in cfg.rng(0,N) :                                                             #for a=0:N
        for b in cfg.rng(0,N) :                                                         #    for b=0:N
            for c in cfg.rng(0,N) :                                                     #        for c=0:N
                temp = 0                                                                #            temp=0;
                for i2 in cfg.rng(0,a) :                                                #            for i2=0:a
                    for j2 in cfg.rng(0,b) :                                            #                for j2=0:b
                        for k2 in cfg.rng(0,c) :                                        #                    for k2=0:
                            temp += C[1,i2,j2,k2]*C[2,a-i2,b-j2,c-k2] #!                #                        temp=temp+C(2,i2+1,j2+1,k2+1)*C(3,a-i2+1,b-j2+1,c-k2+1);
                            #f.write('{} {}\n'.format(C[1,i2,j2,k2],C[2,a-i2,b-j2,c-k2]))
                                                                                        #                    end
                                                                                        #                end
                                                                                        #            end          
                D[a,b,c] = temp #!                                                      #            D(a+1,b+1,c+1)=temp;
                                                                                        #        end
                                                                                        #    end
                                                                                        #end
    return D    

#@profilehooks.timecall
def generate_indices(N) :
    a,b,c,i2,j2,k2 = numpy.mgrid[0:N+1,0:N+1,0:N+1,0:N+1,0:N+1,0:N+1]
    mask = (i2 <= a) & ( j2 <= b ) & ( k2 <= c )
    #return a[mask],b[mask],c[mask],i2[mask],j2[mask],k2[mask]
    return a,b,c,i2,j2,k2,mask

#@profilehooks.timecall
def thesum(D) : return numpy.sum(D,axis=(3,4,5))

def bar(a,b,c,i2,j2,k2) : return a-i2, b-j2, c-k2

def foo(C,a,b,c,i2,j2,k2) :
    #one = numpy.ones_like(a)
    #two = 2*one
    i3,j3,k3 = bar(a,b,c,i2,j2,k2)
    C1, C2 = C[1,:,:,:], C[2,:,:,:]
    return C1[i2,j2,k2], C2[i3,j3,k3]

@profilehooks.profile(filename='vect.prfl')
def vectorized(cfg,C,N) :
    a,b,c,i2,j2,k2,mask = generate_indices(N)
    #C1, C2 = C[1,i2,j2,k2], C[2,a-i2,b-j2,c-k2]
    C1, C2 = foo(C,a,b,c,i2,j2,k2)
    D = C1*C2
    D[~mask] = 0.0
    return thesum(D)

def Dabc_orig(cfg,C,N) :
    vec = vectorized(cfg,C,N)
    foo = itertools_product(cfg,C,N)
    bar = orig(cfg,C,N)
    assert numpy.all(vec == bar)
    sys.exit()
    return bar
