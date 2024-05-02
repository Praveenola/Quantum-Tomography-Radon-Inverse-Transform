%  UFWT_CDJV -- Forward Wavelet Transform (boundary-corrected, NOT preconditioned)
%   Usage
%     wc = UFWT_CDJV(x,L,N)
%   Inputs
%     y    1-d signal, length(x) = 2^J
%     L    Level of V_0;  L << J
%     N    Degree of Daubechies Filters
% 
%   Description
%     CDJV have developed an algorithm for wavelets on the interval which
%     preserves the orthogonality, vanishing moments, smoothness, and compact
%     support of Daubechies wavelets on the line.
% 
%     The algorithm for wavelets on the interval of CDJV involves four objects
%     not present in the usual periodized algorithm: right edge filters, left
%     edge filters, and pre- and post- conditioning operators.
% 
%     These objects are supplied by appropriate requests to MakeCDJVFilter.
% 
% 	 This variant does not apply the preconditioning operator.
% 
%     To reconstruct use CDJV_IWT.
% 
%   See Also
%     IWT_CDJV, FWT_PO, IWT_PO, MakeCDJVFilter
% 

function wc = UFWT_CDJV(x, L, N)
    % Ensure the length of the signal is a power of 2
    if ~isvector(x) || ~isrow(x)
        error('Input signal x must be a 1-D row vector.');
    end

    % Check that L is within a valid range
    if L < 1 || L >= length(x)
        error('Level L is not within a valid range.');
    end

    % Check that N is a valid degree for Daubechies filters
    if N < 1 || N > 20
        error('Degree N is not within a valid range.');
    end

    % Initialize the wavelet coefficients vector
    wc = zeros(size(x));

    % Apply the forward wavelet transform without preconditioning
    wc = CDJV_FWT(x, L, N);
end

function wc = CDJV_FWT(x, L, N)
    % CDJV Forward Wavelet Transform without preconditioning

    J = log2(length(x));

    % Create Daubechies filters
    [h0, h1, g0, g1] = MakeCDJVFilter(N);

    % Initialize the wavelet coefficients
    wc = zeros(size(x));

    % Pre-filtering
    wc = PreFilter(x, h0);

    for j = L:J-1
        % Convolve with high-pass filter
        wc(1:2^(j+1)) = ConvolveHigh(wc(1:2^(j+1)), h1);

        % Upsample
        wc(1:2^(j+1)) = Upsample(wc(1:2^(j+1)), 2);

        % Convolve with low-pass filter
        wc(1:2^(j+1)) = ConvolveLow(wc(1:2^(j+1)), h0);
    end
end

function wc = PreFilter(x, h0)
    % Apply pre-filter to input signal x
    wc = conv(x, h0, 'valid');
end

function y = ConvolveHigh(x, h1)
    % Convolve x with high-pass filter h1
    y = conv(x, h1, 'valid');
end

function y = Upsample(x, L)
    % Upsample signal x by factor L
    y = zeros(1, L * length(x));
    y(1:L:end) = x;
end

function y = ConvolveLow(x, h0)
    % Convolve x with low-pass filter h0
    y = conv(x, h0, 'valid');
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