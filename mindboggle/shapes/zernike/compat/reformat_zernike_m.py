def reformat_zernike(cfg,Z,N) :                                 #function ZM = reformat_zernike(Z,N)
                                                                #
    coord = cfg.load_coord()                                    #load coord;
                                                                #
    ZM = []                                                     #ZM = [];
    for c in cfg.rng(coord,1) :                                 #for c = 1:size(coord,1);
                                                                #    
                                                                #    
        if coord[c,1]+1 <= N+1 :                                #    if coord(c,1)+1 <= N+1;
            ZM[c,1] = Z[coord[c,1]+1,coord[c,2]+1,coord[c,3]+1] #        ZM(c,1) = Z(coord(c,1)+1,coord(c,2)+1,coord(c,3)+1);
        else :                                                  #    else
            break                                               #        break
                                                                #    end
                                                                #end
    return ZM
