%  Inv_PseudopolarFFT: Generalized inverse of Pseudopolar FFT
%   Usage:
%      X = Inv_PseudopolarFFT(Y,Age,Pre,DC,K,E)
%   Inputs:
%     Y      2n*2n matrix (theta,r)
%     Age	Flag (1/0) for use of old code default=0
%     Pre	Preconditioning flag (1/0) default=0
%     MaxIts max number of iterations. Default 5.
%     ErrTol error tolerance. Default 1.e-9.
%   Outputs:
%     X      n*n matrix (x,y)
%   Description:
%     Performs inverse of pseudopolar
%   FFT on 2n*2n matrix by using frft.
%   This is an approximate inverse.
% 
%  This function finds the inverse polar fft
%  by using a conjugate gradient solver on the associated Gram system.
%  The number of iterations is bounded by MaxIts
%  and the residual error by ErrTol.
% 

function X = Inv_PseudopolarFFT(Y, Age, Pre, MaxIts, ErrTol)
    % Check if Age, Pre, MaxIts, and ErrTol are provided and set defaults if not
    if nargin < 2
        Age = 0;  % Default value for the Age flag
    end
    if nargin < 3
        Pre = 0;  % Default value for the Preconditioning flag
    end
    if nargin < 4
        MaxIts = 5;  % Default maximum number of iterations
    end
    if nargin < 5
        ErrTol = 1e-9;  % Default error tolerance
    end

    % Ensure Y is a 2n*2n matrix
    [n, ~] = size(Y);
    if mod(n, 2) ~= 0
        error('Input matrix Y must be of size 2n*2n.');
    end

    % Initialize the output matrix X
    X = zeros(n, n);

    % Define some constants
    theta = linspace(-pi, pi, n);
    r = linspace(0, sqrt(2), n);
    delta_theta = 2 * pi / n;
    delta_r = sqrt(2) / n;

    for i = 1:n
        for j = 1:n
            % Calculate the coordinates in Cartesian space
            x = r(j) * cos(theta(i));
            y = r(j) * sin(theta(i));

            % Perform the inverse pseudopolar FFT using conjugate gradient
            X(i, j) = ConjGrad(Y, x, y, Age, Pre, MaxIts, ErrTol);
        end
    end
end

function x = ConjGrad(Y, x, y, Age, Pre, MaxIts, ErrTol)
    % This is a placeholder for the conjugate gradient solver.
    % You need to implement the conjugate gradient solver
    % or use an existing solver to approximate the inverse.

    % A simple example:
    %x = Y;  % Replace with the actual conjugate gradient solver
    
    % Flatten the 2D input and output matrices for simplicity
    X = X(:);
    y = X;

    r = r(:);
    theta = theta(:);

    % Calculate the matrix-vector product A*y
    function Ay = ApplyMatrixA(Y)
        % Assuming that A is the identity operator, you can modify this part
        Ay = Y;
    end

    % Initialize vectors and variables
    r0 = r - ApplyMatrixA(y);
    p = r0;
    rsold = r0' * r0;

    for i = 1:MaxIts
        % Calculate Ap
        Ap = ApplyMatrixA(p);

        % Calculate the step size (alpha)
        alpha = rsold / (p' * Ap);

        % Update the solution
        y = y + alpha * p;

        % Calculate the new residual
        r = r0 - alpha * Ap;

        % Calculate the squared norm of the new residual
        rsnew = r' * r;

        % Check for convergence based on the squared error
        if sqrt(rsnew) < ErrTol
            break;
        end

        % Calculate the new search direction (beta)
        beta = rsnew / rsold;

        % Update the search direction
        p = r + beta * p;

        % Update the squared norm of the residual
        rsold = rsnew;
    end

    % If you need to reshape the result back to the original 2D shape, you can do so
    y = reshape(y, size(X));
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