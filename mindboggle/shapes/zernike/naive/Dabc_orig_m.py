'''VALID'''

def Dabc_orig(C,N) :                                                               #function D=Dabc_orig(C,N)
                                                                                   #% Pre-computes the values for D
                                                                                   #% and store them in matrix form
                                                                                   #% to be used for Geometric Moments

    D = zeros( (N+1,N+1,N+1),dtype=float ) #Warn: hackery                                      #D = zeros(N+1,N+1,N+1);
    for a in xrange(N+1) :  #Warn: xrange                                          #for a=0:N
        for b in xrange(N+1) : #Warn: xrange                                       #    for b=0:N
            for c in xrange(N+1) : #Warn: xrange                                   #        for c=0:N
                temp = 0.0                                                         #            temp=0;
                for i2 in xrange(a+1) : #Warn: xrange                              #            for i2=0:a
                    for j2 in xrange(b+1) : #Warn: xrange                          #                for j2=0:b
                        for k2 in xrange(c+1) : #Warn: xrange                      #                    for k2=0:c
                            temp += C[1,i2,j2,k2]*C[2,a-i2,b-j2,c-k2] # Warn: 0-1  #          temp=temp+C(2,i2+1,j2+1,k2+1)*C(3,a-i2+1,b-j2+1,c-k2+1);
                                                                                   #                    end
                                                                                   #                end
                                                                                   #            end
                D[a,b,c] = temp #Warn: 0-1                                         #            D(a+1,b+1,c+1)=temp;
                                                                                   #        end
                                                                                   #    end
                                                                                   #end
    return D
