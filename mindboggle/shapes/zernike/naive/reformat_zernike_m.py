'''VALID'''

def reformat_zernike(Z,N,coord_file) :                                #function ZM = reformat_zernike(Z,N)
                                                                      #
    coord = load(coord_file)['coord'] #Warn 0-1                       #load coord;
    count = sum(coord[:,0] <= N)                                                                  #
    ZM = zeros(count,dtype=float)                                     #ZM = [];
    for c in xrange(coord.shape[0]) :                                 #for c = 1:size(coord,1);
                                                                      #    
        if coord[c,0] <= N : #Warn: 0-1                               #    if coord(c,1)+1 <= N+1;
            print c+1, coord[c,0], Z[coord[c,0],coord[c,1],coord[c,2]]
            ZM[c] = Z[coord[c,0],coord[c,1],coord[c,2]] #Warn: 0-1    #        ZM(c,1) = Z(coord(c,1)+1,coord(c,2)+1,coord(c,3)+1);
        else :                                                        #    else
            break                                                     #        break
                                                                      #    end
                                                                      #end
    return ZM
