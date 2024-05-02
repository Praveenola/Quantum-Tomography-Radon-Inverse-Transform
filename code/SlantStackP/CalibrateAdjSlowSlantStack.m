%  TestSlowSlantStack
%
function CalibrateAdjSlowSlantStack()
    % Generate a test image
    n = 256;  % Size of the image
    X = phantom(n);
    
    % Perform Slant Stack Radon transform
    S = SlowSlantStack(X);
    
    % Perform the adjoint Slant Stack Radon transform to reconstruct the image
    X_reconstructed = Adj_SlowSlantStack(S);
    
    % Display the original and reconstructed images
    figure;
    subplot(1, 2, 1);
    imshow(X, []);
    title('Original Image');
    
    subplot(1, 2, 2);
    imshow(X_reconstructed, []);
    title('Reconstructed Image');
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