%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------- Implements Level 3 of the Assignment (PSYCHO) --------------
% Function that implements the inverse AACoder. x is the output signal.
% The output signal is stored with Fs = 40kHz and has 2 channels 
% Where:
% ACCSeq1: The struct that contains the frameType,winType,chl.frameF 
%   (the left channels mdct coefficients), chr.frameF (the right channels mdct coefficients)
% fnameOut: The desirable output file name 
%%
function x = iAACoder3(AACSeq3,fNameOut)

global short_bands long_bands

% Load the Psychoacoustic Model bands
bands = load('TableB219.mat');
long_bands = bands.B219a;
short_bands = bands.B219b;
long_bands(:,1:3) = long_bands(:,1:3) + 1;
short_bands(:,1:3) = short_bands(:,1:3) + 1;

% Calculate huffLTU matrix
huffLUT = loadLUT();

% Initialize empty signal x 
signal = NaN(length(AACSeq3)*1024,2);

% Counter that helps reconstruct the initial signal
counter = 1025;

% For every frame in the struct
for i = 1:length(AACSeq3)
    
    sfc_left = decodeHuff(AACSeq3(i,1).chl.sfc, 12, huffLUT)';
    sfc_right = decodeHuff(AACSeq3(i,1).chr.sfc, 12, huffLUT)';
    
    if AACSeq3(i,1).frameType == "ESH"
        sfc_left = reshape(sfc_left,41,8);
        sfc_right = reshape(sfc_right,41,8);
    end 
    decoded_left = decodeHuff(AACSeq3(i,1).chl.stream, AACSeq3(i,1).chl.codebook, huffLUT)';
    decoded_right = decodeHuff(AACSeq3(i,1).chr.stream, AACSeq3(i,1).chr.codebook, huffLUT)';
    
    frameFout_left_tns = iAACquantizer(decoded_left, sfc_left, AACSeq3(i,1).chl.G, AACSeq3(i,1).frameType);
    frameFout_right_tns = iAACquantizer(decoded_right, sfc_right, AACSeq3(i,1).chr.G, AACSeq3(i,1).frameType);
    
    % First apply the inverse TNS to each frame
    frameFout_left = iTNS(frameFout_left_tns,AACSeq3(i,1).frameType,AACSeq3(i,1).chl.TNScoeffs);
    frameFout_right = iTNS(frameFout_right_tns,AACSeq3(i,1).frameType,AACSeq3(i,1).chr.TNScoeffs);
    
    frameF = [frameFout_left frameFout_right];
    if AACSeq3(i,1).frameType == "ESH"
        frameF(:,1:2:15) = frameFout_left;
        frameF(:,2:2:16) = frameFout_right;
    end 
    
    % Then inverse filterbank to the corresponding frame
    frameT = ifilterbank(frameF,AACSeq3(i,1).frameType,AACSeq3(i,1).winType);

    % If its the first iteration then the singal x is set to be the whole
    % first frame
    if i == 1
        signal(1:2048,:) = frameT;
        continue;
    end
    % If its not the first iteration then the first 1024 elements are added
    % to the previous frame's last 1024 elements  
    signal(counter:counter+1023,:) = signal(counter:counter+1023,:) + frameT(1:1024,:);
    counter = counter + 1024;
    
    % The other half of the frame is assinged to the next 1024 of the
    % singal
    signal(counter:counter+1023,:) = frameT(1025:end,:);
end

% Ignore the first 1024 samples that were zero padded
signal = signal(1025:end,:);
signal(signal < -1 ) = -1;
signal(signal > 1 ) = 1;

% Write the audio file
audiowrite(fNameOut,signal,48000);

% Return the signal values to variable x
if nargout == 1
    x = signal;
end

end

