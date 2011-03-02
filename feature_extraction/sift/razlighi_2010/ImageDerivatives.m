function [Derivative]= ImageDerivatives(Input, DerivationDirection)

%%% This function computes the regular image derivatives along any
%%% direction and returns the derivative values as the same size as the
%%% original image size. An extra index is added with the same value of the
%%% last index to achive this goal therefore the end index of the
%%% derivatives in any directions are zero.



[M N R P] = size(Input);
PaddedInput = Input;


if DerivationDirection == 'x'
    PaddedInput(1,:,:,:) = Input(1,:,:,:);
    PaddedInput(2:M+1,:,:,:) = Input;
    PaddedInput(M+2,:,:,:) = Input(M,:,:,:);
    
    Derivative = (PaddedInput(3:end,:,:,:) - PaddedInput(1:end-2,:,:,:))/2;

elseif DerivationDirection == 'y'
    PaddedInput(:,1,:,:) = Input(:,1,:,:);
    PaddedInput(:,2:N+1,:,:) = Input;
    PaddedInput(:,N+2,:,:) = Input(:,N,:,:);

    Derivative = (PaddedInput(:,3:end,:,:) - PaddedInput(:,1:end-2,:,:))/2;

elseif DerivationDirection == 'z'
    PaddedInput(:,:,1,:) = Input(:,:,1,:);
    PaddedInput(:,:,2:R+1,:) = Input;
    PaddedInput(:,:,R+2,:) = Input(:,:,R,:);

    Derivative = (PaddedInput(:,:,3:end,:) - PaddedInput(:,:,1:end-2,:))/2;

elseif DerivationDirection == 'w'
    PaddedInput(:,:,:,1) = Input(:,:,:,1);
    PaddedInput(:,:,:,2:P+1) = Input;
    PaddedInput(:,:,:,P+2) = Input(:,:,:,P);

    Derivative = (PaddedInput(:,:,:,3:end) - PaddedInput(:,:,:,1:end-2))/2;

end    
 
