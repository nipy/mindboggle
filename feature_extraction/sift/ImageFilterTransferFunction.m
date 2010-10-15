function TransferFunction = ImageFilterTransferFunction(varargin)

if(nargin < 1)
	error('At least one argument is required.') ;
end


for k=1:2:length(varargin)
	switch varargin{k}
    
        case 'Type'
            switch varargin{k+1}
                case 'Ideal LPF'
                    Type = 0;
          
                case 'Butterworth LPF'
                    Type = 1;
          
                case 'Gaussian LPF'
                    Type = 2;
                
                case 'Ideal BPF'
                    Type = 10;
          
                case 'Butterworth BPF'
                    Type = 11;
          
                case 'Gaussian BPF'
                    Type = 12;
                
                otherwise
                    error(['Unknown type: ' varargin{k+1} '.']) ;
            end
            
        case 'Bandwith'
            Bandwith = varargin{k+1} ;
      
        case 'Filter Order'
            FilterOrder = varargin{k+1} ;

        case 'Image Size'
            Temp1 = size(varargin{k+1});
            Temp2 = varargin{k+1};
            switch Temp1(2) 
                case 2
                    M = Temp2(1);
                    N = Temp2(2);
                    P = 1;
                case 3
                    M = Temp2(1);
                    N = Temp2(2);
                    P = Temp2(2);
            end
            
        case 'Center Frequency'
            CenterFrequency = varargin{k+1} ;
    
        otherwise
            error(['Unknown parameter: ' varargin{k} '.']) ;
  end
end




if Type == 0
    for i=1:M
        for j=1:N
            for k=1:P
                Distance = sqrt((i-M/2)^2+(j-N/2)^2+(k-P/2)^2);
                if (Distance <= Bandwith)
                    TransferFunction(i,j,k) = 1;
                else
                    TransferFunction(i,j,k) = 0;
                end
            end
        end
    end
    
end

if Type == 1
    for i=1:M
        for j=1:N
            for k=1:P
                TransferFunction (i, j, k) = ...
                     (1 + (sqrt((i-M/2)^2+(j-N/2)^2+(k-P/2)^2)/Bandwith)^(2*FilterOrder))^(-1);
                
            end
        end
    end
end

if Type == 2
    for i=1:M
        for j=1:N
            for k=1:P
                TransferFunction (i, j, k) = exp(-((i-M/2)^2+(j-N/2)^2+(k-P/2)^2)/(2*Bandwith^2));
                
            end
        end
    end
    
end


if Type == 11
    for i=1:M
        for j=1:N
            for k=1:P
                TransferFunction_L (i, j, k) = ...
                     (1 + (sqrt((i-M/2)^2+(j-N/2)^2+(k-P/2)^2)/(CenterFrequency - Bandwith/2))^(2*FilterOrder))^(-1);
                TransferFunction_H (i, j, k) = ...
                     (1 + (sqrt((i-M/2)^2+(j-N/2)^2+(k-P/2)^2)/(CenterFrequency + Bandwith/2))^(2*FilterOrder))^(-1);
                
            end
        end
    end
    TransferFunction = TransferFunction_H - TransferFunction_L;
end


