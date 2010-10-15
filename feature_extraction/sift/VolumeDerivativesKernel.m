function H = VolumeDerivativesKernel(Dim, Axes)


%%%%%%%%%%%%%%%%%%%%%% Constructing the basic vectors
hS(1) =  1;
hS(2) =  2;
hS(3) =  1;

hD(1) =  1;
hD(2) =  0;
hD(3) = -1;


%%%%%%%%%%%%%%%%%%%%%% Two Dimensional Dreviatives
if Dim == 2
    if Axes == 'X' 
        for i = 1:3
            for j=1:3
                H(i,j) = hD(i)*hS(j);
            end
        end
    elseif Axes == 'Y'
        for i = 1:3
            for j=1:3
                H(i,j) = hD(j)*hS(i);
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%% Three Dimensional Dreviatives
if Dim ==3
    if Axes == 'X' 
        for i = 1:3
            for j=1:3
                for k=1:3
                    H(i,j,k) = hD(i)*hS(j)*hS(k);
                end
            end
        end
    elseif Axes == 'Y'
        for i = 1:3
            for j=1:3
                for k=1:3
                    H(i,j,k) = hD(j)*hS(i)*hS(k);
                end
            end
        end
    elseif Axes == 'Z'
        for i = 1:3
            for j=1:3
                for k=1:3
                    H(i,j,k) = hD(k)*hS(j)*hS(i);
                end
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%% Four Dimensional Dreviatives
if Dim == 4
    if Axes == 'X' 
        for i = 1:3
            for j=1:3
                for k=1:3
                    for l=1:3
                        H(i,j,k,l) = hD(i)*hS(j)*hS(k)*hS(l);
                    end
                end
            end
        end
    elseif Axes == 'Y'
        for i = 1:3
            for j=1:3
                for k=1:3
                    for l=1:3
                        H(i,j,k,l) = hD(j)*hS(i)*hS(k)*hS(l);
                    end
                end
            end
        end
    elseif Axes == 'Z'
        for i = 1:3
            for j=1:3
                for k=1:3
                    for l=1:3
                        H(i,j,k,l) = hD(k)*hS(j)*hS(i)*hS(l);
                    end
                end
            end
        end
     elseif Axes == 'W'
        for i = 1:3
            for j=1:3
                for k=1:3
                    for l=1:3
                        H(i,j,k,l) = hD(l)*hS(j)*hS(k)*hS(i);
                    end
                end
            end
        end
    end
end
    


    