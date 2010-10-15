function [Sigma dSigma] = GaussianSigma(InitialSigma, Nominal_Sigma, S, StartingOctave, LastOctave)
 

for Octave = StartingOctave:(LastOctave+StartingOctave)
    for s = -1:S+1
		
        Sigma(Octave+1-StartingOctave, s+2) = ((InitialSigma * 2^((s/S))));
        %Without downsampleing the Octave needs to be included
        %Sigma(Octave+1-StartingOctave, s+2) = ((InitialSigma * 2^(Octave+(s/S))));
        
        if (s+2 == 1 && Octave+1-StartingOctave ==1)
            PreviousSigma = Nominal_Sigma/(2^StartingOctave);
        elseif (s+2 == 1)
            PreviousSigma = Sigma(Octave+1-StartingOctave, s+2);
        else
            PreviousSigma = Sigma(Octave+1-StartingOctave, s+2-1);
        end

        dSigma(Octave+1-StartingOctave, s+2) = sqrt((Sigma(Octave+1-StartingOctave, s+2))^2 - PreviousSigma^2);
    end

end
