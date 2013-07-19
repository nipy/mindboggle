                                                                                            #%% ZERNIKE 3D FROM GEOMETRIC MOMENTS
                                                                                            #
def zernike(cfg,G,N) :                                                                      #function Z=zernike(G,N)
                                                                                            #
                                                                                            #% Compute the 3D Zernike moments
                                                                                            #% from the geometric moments
    V=cfg.zeros(N+1,N+1,N+1,dtype=complex)                                                  #V=zeros(N+1,N+1,N+1);
    W=cfg.zeros(N+1,N+1,N+1,dtype=complex)                                                  #W=zeros(N+1,N+1,N+1);
    X=cfg.zeros(N+1,N+1,N+1,dtype=complex)                                                  #X=zeros(N+1,N+1,N+1);
    Y=cfg.zeros(N+1,N+1,N+1,dtype=complex)                                                  #Y=zeros(N+1,N+1,N+1);
    Z=cfg.zeros(N+1,N+1,N+1,dtype=complex)                                                  #Z=zeros(N+1,N+1,N+1);
                                                                                            #
                                                                                            #%% FOR a,b,c
                                                                                            #% Computing V
    i=cfg.sqrt(-1);                                                                         #i=sqrt(-1);
    for a in cfg.rng(0,cfg.floor(N/2)) :                                                                    #for a=0:floor(N/2)
        for b in cfg.rng(0,N-2*a) :                                                                   #    for b=0:(N-2*a)
            for c in cfg.rng(0,N-2*a-b) :                                                             #        for c=0:(N-2*a-b)
                                                                                            #            % Calculation for a given a,b,c
                tmp=0                                                                       #            tmp=0;
                for alpha in cfg.rng(0,a+c) :                                                         #            for alpha=0:(a+c)
                    tmp=tmp+cfg.power(i,alpha)*cfg.nchoosek(a+c,alpha)*G[2*a+c-alpha,alpha,b] #!  #                tmp=tmp+power(i,alpha)*nchoosek(a+c,alpha)*G(2*a+c-alpha+1,alpha+1,b+1);
                                                                                            #            end
                    V[a,b,c]=tmp #!                                                      #            V(a+1,b+1,c+1)=tmp;
                                                                                            #        end
                                                                                            #    end
                                                                                            #end
                                                                                            #
                                                                                            #% Computing W
    i=cfg.sqrt(-1)                                                                              #i=sqrt(-1);
    for a in cfg.rng(0,cfg.floor(N/2)) :                                                                    #for a=0:floor(N/2)
        for b in cfg.rng(0,N-2*a) :                                                                   #    for b=0:(N-2*a)
            for c in cfg.rng(0,N-2*a-b) :                                                             #        for c=0:(N-2*a-b)
                                                                                            #            % Calculation for a given a,b,c
                tmp=0                                                                       #            tmp=0;
                for alpha in cfg.rng(0,a) :                                                             #            for alpha=0:a
                    tmp=tmp+cfg.power(-1,alpha)*cfg.power(2,a-alpha)*cfg.nchoosek(a,alpha)*V[a-alpha,b,c+2*alpha] #! #                tmp=tmp+power(-1,alpha)*power(2,a-alpha)*nchoosek(a,alpha)*V(a-alpha+1,b+1,c+2*alpha+1);
                                                                                            #            end
                W[a,b,c]=tmp #!                                                         #            W(a+1,b+1,c+1)=tmp;
                                                                                            #        end
                                                                                            #    end
                                                                                            #end
                                                                                            # 
                                                                                            #% Computing X
    i=cfg.sqrt(-1)                                                                             #i=sqrt(-1);
    for a in cfg.rng(0,cfg.floor(N/2)) :                                                                    #for a=0:floor(N/2)
        for b in cfg.rng(0,N-2*a) :                                                                   #    for b=0:(N-2*a)
            for c in cfg.rng(0,N-2*a-b) :                                                             #        for c=0:(N-2*a-b)
                                                                                            #            % Calculation for a given a,b,c
                tmp=0                                                                       #            tmp=0;
                for alpha in cfg.rng(0,a) :                                                             #            for alpha=0:a
                    tmp=tmp+cfg.nchoosek(a,alpha)*W[a-alpha,b+2*alpha,c] #!                 #                tmp=tmp+nchoosek(a,alpha)*W(a-alpha+1,b+2*alpha+1,c+1);
                                                                                            #            end
                X[a,b,c]=tmp #!                                                         #            X(a+1,b+1,c+1)=tmp;
                                                                                            #        end
                                                                                            #    end
                                                                                            #end
                                                                                            #
                                                                                            #
                                                                                            #%% FOR n,l,m,nu
                                                                                            #% Computing Y
    for l in cfg.rng(0,N) :                                                                             #for l=0:N
        for nu in cfg.rng(0,cfg.floor((N-l)/2)) :                                                           #    for nu=0:floor((N-l)/2)
            for m in cfg.rng(0,l) :                                                                     #        for m=0:l
                                                                                            #            % Calculation for a given l,nu,m
                tmp=0                                                                       #            tmp=0;
                for j in cfg.rng(0,cfg.floor((l-m)/2)) :                                                    #            for j=0:floor((l-m)/2)
                    tmp=tmp+cfg.Yljm(l,j,m)*X[nu+j,l-m-2*j,m] #!                            #                tmp=tmp+Yljm(l,j,m)*X(nu+j+1,l-m-2*j+1,m+1);
                                                                                            #            end
                Y[l,nu,m]=tmp #!                                                        #            Y(l+1,nu+1,m+1)=tmp;
                                                                                            #        end
                                                                                            #    end
                                                                                            #end
                                                                                            #
                                                                                            #%Computing Zernike moments#%Computing Zernike moments
    for n in cfg.rng(0,N) :                                                                             #for n=0:N
        for l in cfg.rng(0,n) :                                                                         #    for l=0:n
            if cfg.mod(n-l,2)==0 :                                                              #        if mod(n-l,2)==0
                for m in cfg.rng(0,l) :                                                                 #            for m=0:l
                                                                                            #                % actual computation
                    tmp=0                                                                   #                tmp=0;
                    k=(n-l)/2                                                               #                k=(n-l)/2;
                    for nu in cfg.rng(0,k) :                                                            #                for nu=0:k
                        tmp=tmp+cfg.Qklnu(k,l,nu)*cfg.conj(Y[l,nu,m]) #!                        #                    tmp=tmp+Qklnu(k,l,nu)*conj(Y(l+1,nu+1,m+1));
                                                                                            #                end
                    Z[n,l,m]=(3/(4*cfg.pi))*tmp #!                                          #                Z(n+1,l+1,m+1)=(3/(4*pi))*tmp;
                                                                                            #                
                                                                                            #            end
                                                                                            #        end
                                                                                            #    end
                                                                                            #end
                                                                                            #
                                                                                            #
                                                                                            #%correction for improper sign problem in sum([n l m]) even terms
    for n in cfg.rng(0,N) :                                                                             #for n=0:N
        for l in cfg.rng(0,n) :                                                                         #    for l=0:n
            for m in cfg.rng(0,l) :                                                                     #        for m=0:l
                                                                                            #            % actual computation
                if cfg.mod(cfg.sum([n,l,m]),2) == 0 :                                               #            if mod(sum([n l m]),2) == 0
                    Z[n,l,m] = cfg.real(Z[n,l,m]) - cfg.imag(Z[n,l,m])*i #!         #                Z(n+1,l+1,m+1) = real(Z(n+1,l+1,m+1)) - imag(Z(n+1,l+1,m+1))*i;
                                                                                            #            end
                                                                                            #            
                if cfg.mod(cfg.sum([n,l,m]),2) == 1 :                                               #            if mod(sum([n l m]),2) == 1
                    Z[n,l,m] = -cfg.real(Z[n,l,m]) + cfg.imag(Z[n,l,m])*i #!         #                Z(n+1,l+1,m+1) = -real(Z(n+1,l+1,m+1)) + imag(Z(n+1,l+1,m+1))*i;
                                                                                            #            end
                                                                                            #            
                                                                                            #        end
                                                                                            #    end
                                                                                            #end
    return Z
