%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------- Implements Level 3 of the Assignment (PSYCHO) --------------
% Function that calculates the SNR from the input signal, the bitrate and 
% the compression achieved. 
% was produced by the AACoder and iAACoder.
% Where:
% fnameIn: The path (or the name if its in the same folder) of the .wav
% file
% fnameOut: The desirable output file name 
%%
function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, fnameAACoded)

% Read the initial signal to calculate the SNR 
[y,fs] = audioread(fNameIn,'double');
y_info = audioinfo(fNameIn);
% Zero pad at the end of the signal
y_padded = [y;zeros(1024,2)];

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
SNR = snr(y_padded,x(1:length(y_padded),:) - y_padded);

% Calculate number of bits for the compressed file 
total_bits = 0;

% Calculation of the total bits of the compressed signal
% The window types require 1 bit, the frameType 2 bits, the codebook 4
% bits,the tns coefficients 4 bitw, the global gain 8 bits, and the stream
% and sfc bits are equal to their lengths as returned from huffman
for i = 1:length(AACSeq3)
    
    window_bits = 1;
    frameType_bits = 2;
    sfc_bits = length(AACSeq3(1).chl.sfc) + length(AACSeq3(1).chr.sfc);
    stream_bits = length(AACSeq3(1).chl.stream) + length(AACSeq3(1).chr.stream);
    tns_coeffs_bits = (length(AACSeq3(1).chl.TNScoeffs) +  length(AACSeq3(1).chr.TNScoeffs)) * 4;
    codebook_bits = (length(AACSeq3(1).chl.codebook) + length(AACSeq3(1).chr.codebook)) * 4;
    global_gain_bits = (length(AACSeq3(1).chl.G) + length(AACSeq3(1).chr.G)) * 8;
    total_bits = total_bits + sfc_bits + stream_bits + tns_coeffs_bits + codebook_bits + global_gain_bits + window_bits + frameType_bits;
    
end

% Compressed bitrate
bitrate = fs * total_bits/(2*length(x));

% Bitrate for the initial signal 
bitrate_initial = y_info.SampleRate * y_info.BitsPerSample * 2;

% Compression
compression = bitrate_initial / bitrate;
end

