% 
%                                                         
% This function computes a fractional Fourier            
%                                  transformation in 1-D 
%                                                        
%

function X_hat = frft(x, alpha)
    % Ensure that the input signal length is a power of 2 (dyadic)
    n = length(x);
    if bitand(n, n - 1) ~= 0
        error('Input signal length must be a power of 2 (dyadic).');
    end
    
    X_hat = fft(x);
    
    % Apply the fractional Fourier transform (FRFT) with angle alpha
    for k = 0:n-1
        X_hat(k+1) = X_hat(k+1) * exp(-1i * alpha * pi * k * (k - n) / (2 * n));
    end
    
    X_hat = ifft(X_hat);
end

%
% Part of BeamLab Version:200
% Built:Friday,23-Aug-2002 00:00:00
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail beamlab@stat.stanford.edu
%
%
% Part of BeamLab Version:200
% Built:Saturday,14-Sep-2002 00:00:00
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail beamlab@stat.stanford.edu
%