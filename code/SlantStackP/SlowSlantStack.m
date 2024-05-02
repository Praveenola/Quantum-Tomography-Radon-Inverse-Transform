%  SlowSlantStack: Slant Stack Radon by Direct Algorithm
%   Usage:
%      S = SlowSlantStack(X)
%   Inputs:
%     X      n*n matrix (x,y)
%   Outputs:
%     S      2n*2n matrix (t,theta)
%   Description:
%     Sums array along lines.
%   See Also:
% 	 FastSlantStack, Adj_FastSlantStack, Inv_FastSlantStack
% 

function S = SlowSlantStack(X)
    % Get the size of the input matrix X
    [n, ~] = size(X);
    
    % Calculate the number of rows and columns for the output matrix S
    m = 2 * n;

    % Initialize the output matrix S
    S = zeros(m, m);

    % Define some constants
    delta_t = 2 * n / m;
    delta_theta = pi / m;

    for i = 1:m
        for j = 1:m
            % Calculate the coordinates in Slant Stack Radon space
            t = (j - 1) * delta_t - n;
            theta = (i - 1) * delta_theta - pi;

            % Sum the input array along the lines in Slant Stack Radon space
            S(i, j) = SlantStackRadonSum(X, t, theta);
        end
    end
end

function sum_val = SlantStackRadonSum(X, t, theta)
    % Calculate the sum along the line defined by (t, theta) in the input array X.
    % This involves traversing the line and summing the values of X along it.

    % Implement the summation algorithm for Slant Stack Radon transform.
    % This will depend on the specific line traversal method you want to use.

    % A simple example:
    %sum_val = sum(X(:));  % Replace with the actual summation along the line

    % Calculate the sum along the line defined by (t, theta) in the input array X.
    
    % Calculate the slope of the line
    slope = tan(theta);
    
    % Initialize the sum to zero
    sum_val = 0;

    % Determine whether the line is more horizontal or vertical
    if abs(slope) < 1
        % Line is more horizontal
        for y = 1:size(X, 2)
            x = t + slope * (y - 1);
            % Check if the coordinates are within the bounds of X
            if x >= 1 && x <= size(X, 1)
                % Linear interpolation for sub-pixel accuracy
                x1 = floor(x);
                x2 = x1 + 1;
                if x2 <= size(X, 1)
                    weight1 = x2 - x;
                    weight2 = x - x1;
                    % Sum the weighted values
                    sum_val = sum_val + weight1 * X(x1, y) + weight2 * X(x2, y);
                end
            end
        end
    else
        % Line is more vertical
        for x = 1:size(X, 1)
            y = (x - t) / slope + 1;
            % Check if the coordinates are within the bounds of X
            if y >= 1 && y <= size(X, 2)
                % Linear interpolation for sub-pixel accuracy
                y1 = floor(y);
                y2 = y1 + 1;
                if y2 <= size(X, 2)
                    weight1 = y2 - y;
                    weight2 = y - y1;
                    % Sum the weighted values
                    sum_val = sum_val + weight1 * X(x, y1) + weight2 * X(x, y2);
                end
            end
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