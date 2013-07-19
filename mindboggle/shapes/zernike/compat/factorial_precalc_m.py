def factorial_precalc(cfg,N) :                         #function F = factorial_precalc(N)
                                                       #
    F = cfg.zeros(N+3)
    for i in cfg.rng(1,N) :                            #for i = 1:N
                                                       #    
        for j in cfg.rng(0,N-i) :                      #    for j=0:(N-i)
            for k in cfg.rng(0,N-i-j) :                #        for k=0:(N-i-j)
                                                       #            
                F[i] = cfg.factorial(i) #!             #            F(i + 1) = factorial(i);
                F[j] = cfg.factorial(j) #!             #            F(j + 1) = factorial(j);
                F[k] = cfg.factorial(k) #!             #            F(k + 1) = factorial(k);
                F[i+j+k+2] = cfg.factorial(i+j+k+2) #! #            F(i+j+k+2 + 1) = factorial(i+j+k+2);
                                                       #            
                                                       #        end
                                                       #    end
                                                       #end
    return F
