%  ifft_mid0: Inverse FFT with grid midpoint at 0
%  Usage:
%    X = ifft_mid0(Y)
%  Inputs:
%    Y	Array(n) 
%  Outputs:
%    X   Array(n)
%  Description:
%   Performs 1d ifft with grid (-n/2):(n/2-1) on both time
%   and frequency side. 
%     y(k) = sum_{t=-n/2}^{n/2-1} exp(i 2pi/n kt) x(t) , (-n/2) <= k < n/2
% 

function X = ifft_mid0(Y)
    % Get the size of the input array Y
    n = length(Y);
    
    % Initialize the output array X
    X = zeros(size(Y));
    
    % Perform the inverse 1D FFT with grid (-n/2):(n/2-1) in both domains
    for t = -n/2:n/2-1
        t_idx = t + n/2 + 1;
        X(t_idx) = 0;
        for k = -n/2:n/2-1
            k_idx = k + n/2 + 1;
            X(t_idx) = X(t_idx) + Y(k_idx) * exp(1i * 2 * pi * k * t / n);
        end
        X(t_idx) = X(t_idx) / n;
    end
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