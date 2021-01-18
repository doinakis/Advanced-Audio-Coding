%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------- Implements Level 2 of the Assignment (TNS) ---------------
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

% Check if the output name icludes the file extension
if length(fNameOut) < 4
    fNameOut = strcat(fNameOut,'.wav');
elseif ~(strcmp(fNameOut(end-3:end),'.wav'))
    fNameOut = strcat(fNameOut,'.wav');
end
% Apply the AAC2 and inverse AAC2
AACSeq2 = AACoder2(fNameIn);
x = iAACoder2(AACSeq2,fNameOut);

% Calculate the SNR of the signal
SNR = snr(y,x(1:length(y),:) - y);

end

