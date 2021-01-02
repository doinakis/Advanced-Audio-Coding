%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that calculates the SNR from the input signal and the one that 
% was produced by the AACoder and iAACoder.
% Where:
% fnameIn: The path (or the name if its in the same folder) of the .wav
% file
% fnameOut: The desirable output file name 
%%
function SNR = demoAAC2(fNameIn, fNameOut)

% Read the initial signal to calculate the SNR 
y = audioread(fNameIn,'double');

% Zero pad at the start of the signal
y = [y;zeros(1024,2)];

% Apply the AAC and inverse AAC 
AACSeq2 = AACoder2(fNameIn);
x = iAACoder2(AACSeq2,fNameOut);

% Calculate the SNR of the signal
SNR = snr(y,x(1:length(y),:) - y);

end

