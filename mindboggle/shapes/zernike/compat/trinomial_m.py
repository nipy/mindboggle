def trinomial(cfg,i,j,k) :                                   #function t=trinomial(i,j,k)
                                                             #% Computes the trinomial of
                                                             #% the three input arguments
                                                             #
    aux_1=cfg.factorial(i+j+k)                               #aux_1=factorial(i+j+k);
    aux_2=cfg.factorial(i)*cfg.factorial(j)*cfg.factorial(k) #aux_2=factorial(i)*factorial(j)*factorial(k);
    t = aux_1/aux_2                                          #t= aux_1/aux_2;
                                                             #
    return t
