%  fftshift1d -- perform 1-d fftshift on columns of 2-d array
%  Usage:
%    Y = fftshift1d(X)
%  Inputs:
%    X	Array(n,m) 
%  Outputs:
%    Y   Array(n,m)
%  Description:
%   Performs fftshift() on each column of X
%

function Y = fftshift1d(X)
    % Get the size of the input array X
    [n, m] = size(X);
    
    % Initialize the output array Y
    Y = zeros(n, m);
    
    % Perform fftshift on each column of X
    for i = 1:m
        Y(:, i) = fftshift(X(:, i));
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