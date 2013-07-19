import profilehooks
import numpy
import numpy.random
import sys
import time
import itertools

def itertools_product(cfg,C,N) :
    D = cfg.zeros(N+1,N+1,N+1)
    for a,b,c in cfg.rng_prod((0,N),repeat=3 ) :
        temp = 0
        for i2,j2,k2 in cfg.rng_prod((0,a),(0,b),(0,c)) :
            temp += C[1,i2,j2,k2]*C[2,a-i2,b-j2,c-k2] #!
        D[a,b,c] = temp #!
    return D

#@profilehooks.timecall(immediate=True)
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
                            #f.write('{} {} {} {} {} {}\n'.format(a,b,c,i2,j2,k2))
                            temp += C[1,i2,j2,k2]*C[2,a-i2,b-j2,c-k2] #!                #                        temp=temp+C(2,i2+1,j2+1,k2+1)*C(3,a-i2+1,b-j2+1,c-k2+1);
                                                                                        #                    end
                                                                                        #                end
                                                                                        #            end          
                D[a,b,c] = temp #!                                                      #            D(a+1,b+1,c+1)=temp;
                                                                                        #        end
                                                                                        #    end
                                                                                        #end
    return D

def generate_indices(N) :
    a,b,c,i2,j2,k2 = numpy.mgrid[0:N+1,0:N+1,0:N+1,0:N+1,0:N+1,0:N+1]
    mask = (i2 <= a) & ( j2 <= b ) & ( k2 <= c )
    return a,b,c,i2,j2,k2,mask

def simple(cfg,C,N) :
    a,b,c,i2,j2,k2,mask = generate_indices(N)
    C1, C2 = C[1,i2,j2,k2], C[2,a-i2,b-j2,c-k2]
    D = C1*C2
    D[~mask] = 0.0
    return numpy.sum(D,axis=(3,4,5))

def nd_reversed(C) : return C[::-1,::-1,::-1]
def iter_rev_sum(cfg,C,N) :
    C1, C2 = C[1,:,:,:], C[2,:,:,:]
    D = numpy.zeros_like(C1)
    for a,b,c in cfg.rng_prod((0,N),repeat=3) :
        c1 = C1[:a+1,:b+1,:c+1]
        c2 = nd_reversed(C2[:a+1,:b+1,:c+1])
        D[a,b,c] = numpy.sum(c1*c2)
    return D

'''
@profilehooks.profile(filename='Dabc_orig.prfl')
def test_profile() :
    from . import MultiprocPipeline as Pipeline
    N = 20
    C = numpy.random.rand(N+1,N+1,N+1,N+1)
    pl = Pipeline()
    #vec = vectorized(pl,C,N)
    #foo = itertools_product(pl,C,N)
    sv = sumversion(pl,C,N)
    o = orig(pl,C,N)
    print numpy.max(numpy.abs(sv-o))
    #assert numpy.all(vec == o)
    #assert numpy.all(foo == o)
    assert numpy.all(sv == o)
'''
