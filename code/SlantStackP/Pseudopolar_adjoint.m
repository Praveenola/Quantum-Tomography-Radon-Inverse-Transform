%  Pseudopolar_adjoint: Adjoint of Pseudopolar FFT
%   Usage:
%      Y=Pseudopolar_adjoint(X)
%   Inputs:
%     X      2n*2n matrix (theta,r)
%   Outputs:
%     Y      n*n matrix (x,y)
%   Description:
%     Performs adjoint of pseudopolar
%   FFT on 2n*2n matrix by using frft.
%   This is an approximate inverse.
%

function Y = Pseudopolar_adjoint(X)
    % Get the size of the input matrix X
    [n, ~] = size(X);
    
    % Ensure that the input matrix X is 2n*2n
    if mod(n, 2) ~= 0
        error('Input matrix X must be of size 2n*2n.');
    end

    % Initialize the output matrix Y
    Y = zeros(n, n);

    % Define some constants
    delta_theta = 2 * pi / n;
    delta_r = sqrt(2) / n;

    for i = 1:n
        for j = 1:n
            % Calculate the coordinates in Cartesian space
            x = -sqrt(2) + (j - 1) * delta_r;
            y = -pi + (i - 1) * delta_theta;

            % Perform the adjoint Pseudopolar FFT using the FRFT
            Y(i, j) = AdjointFRFT(X, x, y);
        end
    end
end

function y = AdjointFRFT(X, x, y)
    % This is a placeholder for the adjoint FRFT operation.
    % You need to implement the adjoint FRFT or use an existing function.

    % A simple example:
    %y = X;  % Replace with the actual adjoint FRFT operation

    % Define the fractional Fourier transform parameter
    alpha = -y;  % Assuming the conjugate of the parameter

    % Perform the adjoint FRFT using the built-in frft function
    Y = frft(X, alpha);

    % Note: In the above code, we assume that the Signal Processing Toolbox
    % is available, which includes the `frft` function for fractional Fourier
    % transforms. If you don't have this toolbox, you can implement the
    % adjoint FRFT using a different method or available library.
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