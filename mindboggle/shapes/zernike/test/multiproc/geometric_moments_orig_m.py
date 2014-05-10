#import profilehooks

#@profilehooks.timecall
def orig(cfg,X,K,N,num_facets,num_vertices) :                                 #function G=geometric_moments_orig(X,K,N,num_facets,num_vertices)
                                                                              #% Computes the geometric moments
                                                                              #% of the volumetric object given
                                                                              #
                                                                              #% PRE-COMPUTATIONS
    cfg.display('Trinomial')                                                  #display('Trinomial')
                                                                              #% Matrix to store the pre-computed trinomial
    tri_matrix=cfg.trinomial_matrix(N)                                        #tri_matrix=trinomial_matrix(N);
                                                                              #% Volume of each tetrahedra
    Vol=cfg.zeros(num_facets,1)                                               #Vol=zeros(num_facets,1);
                                                                              #% Matrix that stored the monomial combinations
    C_temp=cfg.zeros(num_vertices,N+1,N+1,N+1)                                #%C_temp=zeros(num_vertices,N+1,N+1,N+1);
                                                                              #%keyboard
                                                                              #%size(C)
                                                                              #
    cfg.display('Dabc_orig')                                                  #display('Dabc_orig')
                                                                              #%matlabpool close
                                                                              #%matlabpool local 8
    D=cfg.zeros(num_facets,N+1,N+1,N+1)                                       #D=zeros(num_facets,N+1,N+1,N+1);
    C1=cfg.zeros(num_facets,N+1,N+1,N+1)                                      #C1=zeros(num_facets,N+1,N+1,N+1);
    for facet in cfg.rng(1,num_facets) :                                      #for facet=1:num_facets
                                                                              #    % pre-compute D
        cfg.tic()                                                             #    tic
        if cfg.mod(facet,1000) == 0 :                                         #    if mod(facet,1000) == 0
            cfg.display('Facet: {} of {}',facet,num_facets)                   #        display(['Facet: ' num2str(facet) ' of ' num2str(num_facets)])
                                                                              #    end
        C_temp, Vol_temp = cfg.D_CV_orig(facet,num_vertices,N,X,K,tri_matrix) #    [C_temp Vol_temp] = D_CV_orig(facet,num_vertices,N,X,K,tri_matrix);
        D[facet-1,:,:,:] = cfg.Dabc_orig(C_temp,N) #!                         #    D(facet,:,:,:)=Dabc_orig(C_temp,N);
        C1[facet-1,:,:,:] = C_temp[0,:,:,:] #!                                #    C1(facet,:,:,:) = C_temp(1,:,:,:);
        Vol[facet-1,0] = Vol_temp #!                                          #    Vol(facet,1) = Vol_temp;
        if cfg.mod(facet,1000) == 0 :                                         #    if mod(facet,1000) == 0
            cfg.display('Computed facet: {} in {} seconds',facet,cfg.toc())   #        display(['Computed facet: ' num2str(facet) ' in ' num2str(toc) 'seconds'])
                                                                              #    end
                                                                              #end
                                                                              #clear C_temp Vol_temp
                                                                              #
    cfg.display('D_SG_orig')                                                  #display('D_SG_orig')
                                                                              #%matlabpool close
                                                                              #%matlabpool local 8
                                                                              #% GEOMETRIC MOMENTS;
    F = cfg.factorial_precalc(N)                                              #F = factorial_precalc(N);
    G = cfg.zeros(N+1,N+1,N+1)                                                #G=zeros(N+1,N+1,N+1);
    for i in cfg.rng(0,N) :                                                   #for i=0:N
        cfg.tic()                                                             #    tic
        cfg.display('Order: {} of {}',i,N)                                    #    display(['Order: ' num2str(i) ' of ' num2str(N)])
        G[i,:,:] = cfg.D_SG_orig(num_facets,i,N,C1,D,Vol,F) #!                #    G(i+1,:,:)=D_SG_orig(num_facets,i,N,C1,D,Vol,F);
        cfg.display('Computed Order {} in {} seconds',i,cfg.toc())            #    display(['Computed Order ' num2str(i) ' in ' num2str(toc) 'seconds'])
                                                                              #end
    return G

#def D_CV_orig_worker(cfg,facet,num_vertices,N,X,K,tri_matrix) :
#    C_temp,Vol_temp = cfg.D_CV_orig(facet,num_vertices,N,X,K,tri_matrix)
#    return facet,C_temp,Vol_temp

#@profilehooks.timecall
def mproc(cfg,X,K,N,num_facets,num_vertices) :                                 #function G=geometric_moments_orig(X,K,N,num_facets,num_vertices)
    import multiprocessing
                                                                              #% Computes the geometric moments
                                                                              #% of the volumetric object given
                                                                              #
                                                                              #% PRE-COMPUTATIONS
    cfg.display('Trinomial')                                                  #display('Trinomial')
                                                                              #% Matrix to store the pre-computed trinomial
    tri_matrix=cfg.trinomial_matrix(N)                                        #tri_matrix=trinomial_matrix(N);
                                                                              #% Volume of each tetrahedra
    Vol=cfg.zeros(num_facets,1)                                               #Vol=zeros(num_facets,1);
                                                                              #% Matrix that stored the monomial combinations
    C_temp=cfg.zeros(num_vertices,N+1,N+1,N+1)                                #%C_temp=zeros(num_vertices,N+1,N+1,N+1);
                                                                              #%keyboard
                                                                              #%size(C)
                                                                              #
    D=cfg.zeros(num_facets,N+1,N+1,N+1)                                       #D=zeros(num_facets,N+1,N+1,N+1);
    C1=cfg.zeros(num_facets,N+1,N+1,N+1)                                      #C1=zeros(num_facets,N+1,N+1,N+1);
    C_temp_lookup, Vol_temp_lookup = {}, {}
    cfg.display('D_CV_orig')                                                  #display('Dabc_orig')
                                                                              #%matlabpool close
    POOL = multiprocessing.Pool()                                             #%matlabpool local 8
    def D_CV_callback_factory(fct) :
        def callback(result) :
            C_temp, Vol_temp = result
            C_temp_lookup[fct] = C_temp
            Vol_temp_lookup[fct] = Vol_temp
        return callback

    for facet in cfg.rng(1,num_facets) :                                      #for facet=1:num_facets
                                                                                      #    % pre-compute D
        #cfg.tic()                                                             #    tic
        #if cfg.mod(facet,1000) == 0 :                                         #    if mod(facet,1000) == 0
        #    cfg.display('Facet: {} of {}',facet,num_facets)                   #        display(['Facet: ' num2str(facet) ' of ' num2str(num_facets)])
                                                                              #    end
        POOL.apply_async( cfg.D_CV_orig, (facet,num_vertices,N,X,K,tri_matrix), callback=D_CV_callback_factory(facet) ) #    [C_temp Vol_temp] = D_CV_orig(facet,num_vertices,N,X,K,tri_matrix);
    POOL.close()
    POOL.join()

    for facet in cfg.rng(1,num_facets) :
        C1[facet-1,:,:,:] = C_temp_lookup[facet][0,:,:,:] #!                                #    C1(facet,:,:,:) = C_temp(1,:,:,:);
        Vol[facet-1,0] = Vol_temp_lookup[facet] #!                                          #    Vol(facet,1) = Vol_temp;
        #if cfg.mod(facet,1000) == 0 :                                         #    if mod(facet,1000) == 0
            #cfg.display('Computed facet: {} in {} seconds',facet,cfg.toc())   #        display(['Computed facet: ' num2str(facet) ' in ' num2str(toc) 'seconds'])
                                                                              #    end
                                                                              #end
    cfg.display('Dabc_orig')
    POOL = multiprocessing.Pool()
    def Dabc_callback_factory(fct) :
        def callback(result) :
            D[fct-1,:,:,:] = result #!
        return callback

    for facet in cfg.rng(1,num_facets) :
        POOL.apply_async( cfg.Dabc_orig, (C_temp_lookup[facet], N), callback=Dabc_callback_factory(facet))        #    D(facet,:,:,:)=Dabc_orig(C_temp,N);              
    POOL.close()
    POOL.join()
                                                                              #clear C_temp Vol_temp
                                                                              #
    cfg.display('D_SG_orig')                                                  #display('D_SG_orig')
                                                                              #%matlabpool close
                                                                              #%matlabpool local 8
                                                                              #% GEOMETRIC MOMENTS;
    F = cfg.factorial_precalc(N)                                              #F = factorial_precalc(N);
    G = cfg.zeros(N+1,N+1,N+1)                                                #G=zeros(N+1,N+1,N+1);
    POOL = multiprocessing.Pool()
    def D_SG_callback_factory(ii) :
        def callback(result) :
            G[ii,:,:] = result
        return callback

    for i in cfg.rng(0,N) :                                                   #for i=0:N
        POOL.apply_async( cfg.D_SG_orig, (num_facets,i,N,C1,D,Vol,F), callback=D_SG_callback_factory(i)) #!                #    G(i+1,:,:)=D_SG_orig(num_facets,i,N,C1,D,Vol,F);
    POOL.close()
    POOL.join()
                                                                              #end
    return G


#@profilehooks.timecall
def geometric_moments_orig(cfg,X,K,N,num_facets,num_vertices) :               #function G=geometric_moments_orig(X,K,N,num_facets,num_vertices)
    import numpy

    bar = mproc(cfg,X,K,N,num_facets,num_vertices)
    foo = orig(cfg,X,K,N,num_facets,num_vertices)
    assert numpy.all(foo==bar)
    return foo
