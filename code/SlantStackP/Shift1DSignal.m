%  Shift1DSignal -- Circulant shift with fractional shifts allowed
% 	Usage
% 		y = Shift1DSignal(x,delta)
% 	Inputs
% 		x	array(n,1) or array(1,n), n dyadic
% 	Outputs
% 		y	array shaped like x
%

function y = Shift1DSignal(x, delta)
    % Ensure that delta is within the range [0, 1)
    delta = mod(delta, 1);

    % Get the length of the input signal x
    n = length(x);

    % Determine the integer and fractional parts of the shift
    delta_int = floor(delta);
    delta_frac = delta - delta_int;

    % Circularly shift the signal by the integer part of delta
    if delta_int > 0
        y = circshift(x, [delta_int, 0]);
    elseif delta_int < 0
        y = circshift(x, [n + delta_int, 0]);
    else
        y = x;  % No integer shift
    end

    % Apply fractional shift using linear interpolation
    if delta_frac > 0
        y = (1 - delta_frac) * y + delta_frac * circshift(x, [1, 0]);
    elseif delta_frac < 0
        y = (1 + delta_frac) * y - delta_frac * circshift(x, [-1, 0]);
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