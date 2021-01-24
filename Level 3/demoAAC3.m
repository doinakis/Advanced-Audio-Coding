%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------- Implements Level 3 of the Assignment (PSYCHO) --------------
% Function that calculates the SNR from the input signal and the one that 
% was produced by the AACoder and iAACoder.
% Where:
% fnameIn: The path (or the name if its in the same folder) of the .wav
% file
% fnameOut: The desirable output file name 
%%
function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, fnameAACoded)

% Read the initial signal to calculate the SNR 
y = audioread(fNameIn,'double');

% Zero pad at the end of the signal
y = [y;zeros(1024,2)];

% Check if the output name icludes the file extension
if length(fNameOut) < 4
    fNameOut = strcat(fNameOut,'.wav');
elseif ~(strcmp(fNameOut(end-3:end),'.wav'))
    fNameOut = strcat(fNameOut,'.wav');
end

% Check if the output name icludes the file extension
if length(fnameAACoded) < 4
    fnameAACoded = strcat(fnameAACoded,'.mat');
elseif ~(strcmp(fnameAACoded(end-3:end),'.mat'))
    fnameAACoded = strcat(fnameAACoded,'.mat');
end

% Apply the AAC3 and inverse AAC3
AACSeq3 = AACoder3(fNameIn, fnameAACoded);
x = iAACoder3(AACSeq3,fNameOut);
% Calculate the SNR of the signal
SNR = snr(y,x(1:length(y),:) - y);

end

