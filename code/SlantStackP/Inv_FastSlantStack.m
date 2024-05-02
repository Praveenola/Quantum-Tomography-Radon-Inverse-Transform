%  Inv_FastSlantStack: Inverse Radon Transform for Discrete Data
%   Usage:
%      X = Inv_FastSlantStack(S)
%   Inputs:
%     S      2n*2n matrix (t,theta)
%   Outputs:
%     X      n*n matrix (x,y)
%   Description:
%     Backprojects array along lines.
% 
%

function X = Inv_FastSlantStack(S)
    % Get the size of the input matrix S
    [n, ~] = size(S);
    
    % Initialize the output matrix X
    X = zeros(n, n);
    
    % Define the range of theta (angles) and t (offsets)
    theta = linspace(-pi, pi, n);
    t = linspace(-sqrt(2) * n, sqrt(2) * n, n);
    
    for i = 1:n
        for j = 1:n
            % Calculate the line equation parameters
            a = cos(theta(i));
            b = sin(theta(i));
            c = t(j);
            
            % Calculate the coordinates along the line
            x = linspace(-n, n, n);
            y = (-a * x - c) / b;
            
            % Interpolate the values using griddata and store in X
            X(i, j) = sum(interp2(S, x, y, 'linear', 0));
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