function f = GaussianKernel(Sig, Dim)
 


%Create a filter, h, that can be used to approximate linear camera motion.
%h = fspecial('gaussian', FilterSize, Sig);

%iter = (FilterSize-1)/2; 
%iter = ceil(3*Sig);

if Dim < 3
     iter = ceil(4*Sig);   
	 for x = -iter : iter
		for y = -iter : iter
			f(x+iter+1, y+iter+1) = ((2*pi*Sig^2)^-1)*exp(-(x^2 + y^2 )/(2*Sig^2));
		end
	end
elseif Dim > 3
	iter = ceil(3*Sig);
    for x = -iter : iter
		for y = -iter : iter
			for z = -iter : iter
				f(x+iter+1, y+iter+1, z+iter+1) = ((2*pi*Sig^2)^-1.5)*exp(-(x^2 + y^2 + z^2)/(2*Sig^2));
			end
        end
    end
    
end

f = f/sum(f(:));





    