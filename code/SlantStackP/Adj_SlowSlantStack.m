%  Adj_SlowSlantStack: Adjoint Slant Stack Radon by Direct Algorithm
%   Usage:
%      X = Adj_SlowSlantStack(S)
%   Inputs:
%     S      2n*2n matrix (t,theta)
%   Outputs:
%     X      n*n matrix (x,y)
%   Description:
%     Sums array along lines.
%   See Also:
% 	 SlowSlantStack, FastSlantStack, Adj_FastSlantStack, Inv_FastSlantStack
% 
%

function X = Adj_SlowSlantStack(S)
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
            
            % Perform backprojection by summing values along the line
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