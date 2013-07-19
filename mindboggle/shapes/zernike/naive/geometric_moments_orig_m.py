'''UNVALIDATED'''

#from numpy import zeros, mod, savez
#from scipy.io import loadmat
#import time

def geometric_moments_orig(X,K,N,num_facets,num_vertices) :                             #function G = geometric_moments_orig(X,K,N,num_facets,num_vertices)
                                                                                        #% Computes the geometric moments
                                                                                        #% of the volumetric object given
    
                                                                                        #% PRE-COMPUTATIONS
    print 'Trinomial'                                                                   #display('Trinomial')
                                                                                        #% Matrix to store the pre-computed trinomial
    tri_matrix = trinomial_matrix(N)                                                    #tri_matrix = trinomial_matrix(N);
                                                                                        #%Volume of each tetrahedra
    Vol = zeros( (num_facets,1) )                                                       #Vol = zeros(num_facets,1);
                                                                                        #% Matrix that stored the monomial combinations
                                                                                        #%C_temp = zeros(num_vertices,N+1,N+1,N+1);
                                                                                        #%keyboard
                                                                                        #%size(C)

    print 'Dabc_orig'                                                                   #display('Dabc_orig')
                                                                                        #%matlabpool close
                                                                                        #%matlabpool local 8
    D = zeros( (num_facets,N+1,N+1,N+1) ) #Warn: 0-1                                          #D=zeros(num_facets,N+1,N+1,N+1);
    C1 = zeros( (num_facets,N+1,N+1,N+1) ) #Warn: 0-1                                         #C1=zeros(num_facets,N+1,N+1,N+1);
    for facet in xrange(num_facets) :                                                   #for facet=1:num_facets
        print facet,num_facets                                                          #% pre-compute D
        tic = time.time()                                                               #tic
        if mod(facet,1000) == 0 :                                                       #if mod(facet,1000) == 0
            print 'Facet: {} of {}'.format(facet,num_facets)                            #display(['Facet: ' num2str(facet) ' of ' num2str(num_facets)])
                                                                                        #end
        C_temp, Vol_temp = D_CV_orig(facet,num_vertices,N,X,K,tri_matrix)               #[C_temp Vol_temp] = D_CV_orig(facet,num_vertices,N,X,K,tri_matrix);
        D[facet,:,:,:] = Dabc_orig(C_temp,N)                                            #D(facet,:,:,:)=Dabc_orig(C_temp,N);
        C1[facet,:,:,:] = C_temp[0,:,:,:] #Warn: 0-1                                    #C1(facet,:,:,:) = C_temp(1,:,:,:);
        Vol[facet,0] = Vol_temp #Warn: 0-1                                              #Vol(facet,1) = Vol_temp;
        if mod(facet,1000) == 0 :                                                       #if mod(facet,1000) == 0
            toc = time.time() - tic
            print 'Computed facet: {} in {}'.format(facet,toc)                          #display(['Computed facet: ' num2str(facet) ' in ' num2str(toc) 'seconds'])
                                                                                        #end
                                                                                        #end
    del C_temp, Vol_temp                                                                #clear C_temp Vol_temp
    
    print 'D_SG_orig'                                                                   #display('D_SG_orig')
                                                                                        #%matlabpool close
                                                                                        #%matlabpool local 8
    #savez('save_file.npy', D=D,C1=C1,Vol=Vol )
    #raw_input()
    #f = loadmat('../matlab/demo/save_file.mat')
    #D, C1, Vol = f['D'], f['C1'], f['Vol']
                                                                                        #% GEOMETRIC MOMENTS;
    F = factorial_precalc(N)                                                            #F = factorial_precalc(N);
    G = zeros((N+1,N+1,N+1))                                                            #G=zeros(N+1,N+1,N+1);
    for i in xrange(N+1) :                                                                #for i=0:N
        tic = time.time()                                                               #tic
        print 'Order: {} of {}'.format(i,N)                                             #display(['Order: ' num2str(i) ' of ' num2str(N)])
        G[i,:,:] = D_SG_orig(num_facets,i,N,C1,D,Vol,F) #Warn: 0-1                      #G(i+1,:,:)=D_SG_orig(num_facets,i,N,C1,D,Vol,F);
        toc = time.time() - tic
        print 'Computed Order {} in {}'.format(i,toc)                                   #display(['Computed Order ' num2str(i) ' in ' num2str(toc) 'seconds'])
                                                                                        #end
    return G
