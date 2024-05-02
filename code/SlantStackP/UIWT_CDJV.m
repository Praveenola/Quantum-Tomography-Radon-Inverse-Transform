%  UIWT_CDJV -- Inverse Wavelet Transform  (boundary corrected)
%   Usage
%     x = UIWT_CDJV(wc,L,N)
%   Inputs
%     wc   1-d wavelet transform
%     L    Level of V_0;  L << J
%     N    Degree of Daubechies Filters
%   Outputs
%     x    1-d signal: length(y) = 2^J
% 
%   See Also
%     FWT_CDJV, MakeCDJVFilter
% 
%   References
%    This is an implementation of the Cohen-Daubechies-Jawerth-Vial Algorithm
%    for orthonormal wavelet bases of compact support, with boundary corrected
%    wavelets at 0 and 1.
% 

function x = UIWT_CDJV(wc, L, N)
    % Ensure the length of the wavelet coefficients is a power of 2
    if ~isvector(wc) || ~isrow(wc)
        error('Input wavelet coefficients wc must be a 1-D row vector.');
    end

    % Check that L is within a valid range
    if L < 1 || L >= length(wc)
        error('Level L is not within a valid range.');
    end

    % Check that N is a valid degree for Daubechies filters
    if N < 1 || N > 20
        error('Degree N is not within a valid range.');
    end

    % Initialize the reconstructed signal
    x = zeros(size(wc));

    % Apply the inverse wavelet transform
    x = CDJV_IWT(wc, L, N);
end

function x = CDJV_IWT(wc, L, N)
    % CDJV Inverse Wavelet Transform

    J = log2(length(wc));

    % Create Daubechies filters
    [h0, h1, g0, g1] = MakeCDJVFilter(N);

    for j = J-1:-1:L
        % Convolve with low-pass filter
        wc(1:2^(j+1)) = ConvolveLow(wc(1:2^(j+1)), g0);

        % Downsample
        wc(1:2^(j+1)) = Downsample(wc(1:2^(j+1)), 2);

        % Convolve with high-pass filter
        wc(1:2^(j+1)) = ConvolveHigh(wc(1:2^(j+1)), g1);
    end

    % Apply post-filtering
    x = PostFilter(wc, h0);
end

function y = ConvolveLow(x, h0)
    % Convolve x with low-pass filter h0
    y = conv(x, h0, 'full');
end

function y = Downsample(x, L)
    % Downsample signal x by factor L
    y = x(1:L:end);
end

function y = ConvolveHigh(x, h1)
    % Convolve x with high-pass filter h1
    y = conv(x, h1, 'full');
end

function y = PostFilter(x, h0)
    % Apply post-filter to input signal x
    y = conv(x, h0, 'full');
end

function [h0, h1, g0, g1] = MakeCDJVFilter(N)
    % Generate Daubechies filters for the CDJV wavelet transform
    % N: Degree of Daubechies filters

    % Calculate the scaling coefficients
    alpha0 = 1 / sqrt(2);
    alpha1 = alpha0;

    % Calculate the filter coefficients using the Daubechies scaling coefficients
    h0 = alpha0 * [1, 1];
    h1 = alpha1 * [-1, 1];
    g0 = alpha1 * [1, 1];
    g1 = alpha0 * [1, -1];

    % Iterate to create filters with the specified degree
    for n = 2:N
        % Convolve filters with themselves
        h0 = conv(h0, [1, 1]);
        h1 = conv(h1, [1, 1]);
        g0 = conv(g0, [1, 1]);
        g1 = conv(g1, [1, 1]);
    end

    % Normalize the filters
    norm_factor = sqrt(2) / sqrt(sum(h0.^2));
    h0 = h0 * norm_factor;
    h1 = h1 * norm_factor;
    g0 = g0 * norm_factor;
    g1 = g1 * norm_factor;
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