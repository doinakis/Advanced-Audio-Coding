%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------- Implements Level 1 of the Assignment (SSC & filterbank) --------
% Function that implements the inverse AACoder. x is the output signal.
% The output signal is stored with Fs = 40kHz and has 2 channels
% Where:
% ACCSeq1: The struct that contains the frameType,winType,chl.frameF
%   (the left channels mdct coefficients), chr.frameF (the right channels mdct coefficients)
% fnameOut: The desirable output file name
%%
function x = iAACoder1(AACSeq1,fNameOut)

% Initialize empty signal x
signal = NaN(length(AACSeq1)*1024,2);

% Counter that helps reconstruct the initial signal
counter = 1025;

% For every frame in the struct
for i = 1:length(AACSeq1)

    % Initialize a frame with 2 channels
    frameF = [AACSeq1(i,1).chl.frameF AACSeq1(i,1).chr.frameF];

    if AACSeq1(i,1).frameType == "ESH"
        % If the frame is ESH make the frame 128-by-16
        frameF(:,1:2:15) = AACSeq1(i,1).chl.frameF;
        frameF(:,2:2:16) = AACSeq1(i,1).chr.frameF;
    end
    % Apply the inverse filterbank to the corresponding frame
    frameT = ifilterbank(frameF,AACSeq1(i,1).frameType,AACSeq1(i,1).winType);

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

% Write the audio file
audiowrite(fNameOut,signal,48000);

% Return the signal values to variable x
if nargout == 1
    x = signal;
end

end

