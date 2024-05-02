%  FractionalSlowFT: Fractional Fourier Transform (midpoint at 0)
%   Usage:
%     xhat = FractionalSlowFT(x,alpha)
%   Inputs:
%      x       array(n), n dyadic
%      alpha   scalar -- fractional exponent   
%   Outputs:
%      xhat    array(n), n dyadic
%   Description:
%      Performs FT of series x with an nonstandard "alpha" factor
%   in exponent, so matrix element W_{j,k} = e^(-alpha*i*2*pi*k*t/n).
%       alpha=1   corresponds to ordinary FT.
%       alpha=-1  corresponds to (unnormalized) ordinary inverse FT.
%   Transforms with +/-alpha are adjoints of each other
%   Both the inputs and outputs are centered so that
%   k and t run through -n/2 <= k,t < n/2
%

function xhat = FractionalSlowFT(x, alpha)
    % Ensure that the input signal length is a power of 2 (dyadic)
    n = length(x);
    if bitand(n, n - 1) ~= 0
        error('Input signal length must be a power of 2 (dyadic).');
    end
    
    % Initialize the output signal
    xhat = zeros(size(x));
    
    % Calculate the FRFT using a loop
    for t = -n/2:n/2-1
        t_idx = t + n/2 + 1;
        xhat(t_idx) = 0;
        for k = -n/2:n/2-1
            k_idx = k + n/2 + 1;
            xhat(t_idx) = xhat(t_idx) + x(k_idx) * exp(1i * alpha * pi * k * t / n);
        end
    end
    
    % Normalize the result
    xhat = xhat / sqrt(n);
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