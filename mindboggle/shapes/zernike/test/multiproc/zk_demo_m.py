def zk_demo(cfg,V,F,ZMvtk) :                                      #function max_abs_diff = demo( V,F,ZMvtk )
                                                                  #
    num_vertices = 3                                              #  num_vertices = 3;
    num_facets = cfg.size(F,0) #!                                 #  num_facets   = size(F,1);
    N = 20                                                        #  N = 20; %code will only work upto order 20 (can be fixed later)
                                                                  #
                                                                  #  % ZERNIKE MOMENTS
    G = cfg.geometric_moments_orig(V,F,N,num_facets,num_vertices) #  G=geometric_moments_orig(V,F,N,num_facets,num_vertices);
    Z = cfg.zernike(G,N)                                          #  Z=zernike(G,N);
    Descriptors = cfg.feature_extraction(Z,N)                     #  Descriptors=feature_extraction(Z,N);
                                                                  #  ZM= reformat_zernike(Z,N);
                                                                  #
                                                                  #  %compare to your code
                                                                  #  max_abs_diff = max(abs(abs(ZM) - abs(ZMvtk))); %error should be around 1e-8
    return Descriptors
