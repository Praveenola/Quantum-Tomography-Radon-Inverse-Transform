%  FastSlantStack: Radon Transform for Discrete Data
%   Usage:
%      S = FastSlantStack(X)
%   Inputs:
%     X      n*n matrix (x,y)
%   Outputs:
%     S      2n*2n matrix (t,theta)
%   Description:
%     Sums array along lines.
% 

function S = FastSlantStack(X)
    % Get the size of the input matrix X
    [n, ~] = size(X);
    
    % Define the range of theta (angles) and t (offsets)
    theta = linspace(-pi, pi, 2 * n);
    t = linspace(-sqrt(2) * n, sqrt(2) * n, 2 * n);
    
    % Initialize the output matrix S
    S = zeros(2 * n, 2 * n);
    
    for i = 1:2 * n
        for j = 1:2 * n
            % Calculate the line equation parameters
            a = cos(theta(i));
            b = sin(theta(i));
            c = t(j);
            
            % Calculate the values along the line using 2D interpolation
            x = linspace(-n, n, 2 * n);
            y = (-a * x - c) / b;
            
            % Interpolate the values using griddata and store in S
            S(i, j) = sum(interp2(X, x, y, 'linear', 0));
        end
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