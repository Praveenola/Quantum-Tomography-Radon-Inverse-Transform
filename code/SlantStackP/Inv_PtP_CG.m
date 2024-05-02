%  Inv_PtP_CG: Generalized inverse of Pseudopolar FFT Frame
%   Usage:
%      Y = Inv_PtP_CG(X,Age,Pre,MaxIts,ErrTol)
%   Inputs:
%     X      n*n matrix (theta,r)
%     Age	Flag (1/0) for use of old code default=0
%     Pre	Preconditioning flag (2/1/0) default=1
%     MaxIts max number of iterations, default 5
%     ErrTol error tolerance, default 1.e-9
%   Outputs:
%     Y      n*n matrix (x,y)
%   Description:
%     Performs approximate inverse of P'P using conjugate
%     gradients. Usually requires 3 or 4 steps only.
%   See Also
% 	PtP, PseudopolarFFT, Adj_PseudopolarFFT
% 

function Y = Inv_PtP_CG(X, Age, Pre, MaxIts, ErrTol)
    % Check if Age, Pre, MaxIts, and ErrTol are provided and set defaults if not
    if nargin < 2
        Age = 0;    % Default value for the Age flag
    end
    if nargin < 3
        Pre = 1;    % Default value for Preconditioning flag
    end
    if nargin < 4
        MaxIts = 5; % Default maximum number of iterations
    end
    if nargin < 5
        ErrTol = 1e-9;  % Default error tolerance
    end

    % Ensure X is an n*n matrix
    [n, ~] = size(X);
    if mod(n, 2) ~= 0
        error('Input matrix X must be of size n*n.');
    end

    % Initialize the output matrix Y
    Y = zeros(n, n);

    % Define some constants
    delta_theta = 2 * pi / n;
    delta_r = sqrt(2) / n;

    for i = 1:n
        for j = 1:n
            % Calculate the coordinates in Pseudo-Polar space
            theta = -pi + (i - 1) * delta_theta;
            r = j * delta_r;

            % Perform the inverse P'P using conjugate gradient
            Y(i, j) = ConjGrad(X, theta, r, Age, Pre, MaxIts, ErrTol);
        end
    end
end

function y = ConjGrad(X, theta, r, Age, Pre, MaxIts, ErrTol)
    % This is a placeholder for the conjugate gradient solver.
    % You need to implement the conjugate gradient solver
    % or use an existing solver to approximate the inverse.

    % A simple example:
    %y = X;  % Replace with the actual conjugate gradient solver
    
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