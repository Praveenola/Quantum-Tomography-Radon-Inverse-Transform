%  Adj_PseudopolarFFT: Adjoint of Pseudopolar FFT
%   Usage:
%      X = Adj_PseudopolarFFT(Y)
%   Inputs:
%     Y      2n*2n matrix (theta,r)
%     Pre	Preconditioning flag (1/0) default=0
%     DC     Flag for Special Treatment of DC term (1/0) default=0
%   Outputs:
%     X      n*n matrix (x,y)
%   Description:
%     Performs adjoint of pseudopolar FFT.
%   This is an approximate inverse to PseudopolarFFT.
%   See Also
%     PseudopolarFFT, Inv_PseudopolarFFT
%
%

function X = Adj_PseudopolarFFT(Y, Pre, DC)
    % Check if Pre and DC flags are provided, and set defaults if not
    if nargin < 2
        Pre = 0; % Default value for Preconditioning flag
    end
    if nargin < 3
        DC = 0;  % Default value for DC flag
    end
    
    % Get the size of the input matrix Y
    [n, ~] = size(Y);
    
    % Initialize the output matrix X
    X = zeros(n, n);
    
    % Define the range of theta and r values
    theta = linspace(-pi, pi, n);
    r = linspace(0, sqrt(2), n);
    
    % Define constants
    delta_theta = 2 * pi / n;
    delta_r = sqrt(2) / n;
    
    for i = 1:n
        for j = 1:n
            % Calculate the coordinates in Cartesian space
            x = r(j) * cos(theta(i));
            y = r(j) * sin(theta(i));
            
            % Calculate the adjoint of pseudopolar FFT
            if DC == 1 && i == 1 && j == 1
                X(i, j) = sum(sum(Y)) / n^2; % Special treatment for DC term
            else
                X(i, j) = sum(sum(Y .* exp(1i * (x * delta_r + y * delta_theta))) / n^2);
            end
            
            % Apply preconditioning if specified
            if Pre == 1
                X(i, j) = X(i, j) / (r(j) + eps); % Add epsilon to avoid division by zero
            end
        end
    end
end

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