%  PrecondPPFT: Preconditioner for Pseudopolar FFT
%   Usage:
%    PrPPFT = PrecondPPFT(PPFT,Pre,Pow)
%   Inputs:
%    PPFT		array(m,m) PseudopolarFFT
%    Pre			2/1 Preconditioner type
%    Pow         exponent of preconditioner e.g. 1,1/2,-1/2,-1
%   Outputs:
%    PrPPFT      array(m,m) Preconditioned PseudopolarFFT
%   Description:
%    The Pseudopolar FFT array is multiplied by a function of
%     pseudo-radius which is a chosen power of pseudo radius.
%   See Also
%    Pseudopolar FFT, Adj_Pseudopolar FFT
% 

function PrPPFT = PrecondPPFT(PPFT, Pre, Pow)
    % Check the preconditioner type (Pre)
    if Pre == 2
        % Apply preconditioner using the power of pseudo radius
        PrPPFT = PPFT .* (abs(PPFT) .^ Pow);
    elseif Pre == 1
        % Apply a different preconditioner (modify as needed)
        % Example: PrPPFT = PPFT + Pow * randn(size(PPFT));
    else
        error('Invalid preconditioner type. Use 1 or 2.');
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