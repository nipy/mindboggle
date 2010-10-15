function FlipedPoint = PointFliper(InputPoint, FlipDirection, FlipCenter)
 
[M N R] = size(InputPoint);

for i = 1:M
    
    TransformedPoint = InputPoint(i,1:N) - FlipCenter;
    if FlipDirection == 1               %horizental flip
        FlipedTransformedPoint = [-TransformedPoint(1) TransformedPoint(2)]; 
    elseif FlipDirection == 2           %Vertical flip
        FlipedTransformedPoint = [TransformedPoint(1) -TransformedPoint(2)]; 
    end

    FlipedPoint(i,1:N) = FlipedTransformedPoint + FlipCenter;
end