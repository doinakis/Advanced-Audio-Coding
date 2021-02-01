%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------- Implements Level 2 of the Assignment (TNS) ---------------
% Function that calculates the coded frames and returns a K-by-1 srtuct
% with the frameType,winType,chl.TNScoeffs,chr.TNScoeffs (the quantized TNS
% coefficients of the left and right channel, 4-by-8 for ESH 4-by-1
% else), chl.frameF, chr.frameF (the MDCT coefficients after the TNS is
% applied for the left and right channel).
% Where:
% fNameIn: The path (or the name if its in the same folder) of the .wav
% file
%%
function AACSeq2 = AACoder2(fNameIn)

% Load the Psychoacoustic Model bands
global long_bands short_bands
bands = load('TableB219.mat');
long_bands = bands.B219a;
short_bands = bands.B219b;
long_bands(:,1:3) = long_bands(:,1:3) + 1;
short_bands(:,1:3) = short_bands(:,1:3) + 1;

% Read the audio signal from the input
y = audioread(fNameIn,'double');

% Define the window type
window_type = "SIN";

% Preallocate the space for the frames
AACSeq2(ceil(size(y,1)/1024),1) = struct();

% Helping variables to create the 50% overlapping
% one for each channel
frame1 = NaN(2048,ceil(size(y,1)/1024));
frame2 = NaN(2048,ceil(size(y,1)/1024));

% Zero pad at the start and the end of the signal
y = [zeros(1024,2);y;zeros(1024,2)];
frame_counter = 1;

% Create the overlapping frames
for i = 1:1024:size(y,1)

    % If the samples cant completely fill the las frame zero pad at the end
    if i+2047 > size(y,1)
        frame1(1:2048,frame_counter) = [y(i:end,1);zeros(i+2047-size(y,1),1)];
        frame2(1:2048,frame_counter) = [y(i:end,2);zeros(i+2047-size(y,1),1)];
        break;
    end
    frame1(1:2048,frame_counter) = y(i:i+2047,1);
    frame2(1:2048,frame_counter) = y(i:i+2047,2);
    frame_counter = frame_counter + 1;
end


% For every frame apply the AAC coding process
for i = 1:frame_counter

    % Make a dual channel frame from frame1 and frame2
    frameT = [frame1(:,i) frame2(:,i)];

    % For the last frame we assume that the next frame is zeros
    if i == frame_counter
        nextframeT = zeros(2048,2);
    else
        nextframeT = [frame1(:,i+1) frame2(:,i+1)];
    end

    % For the first frame assume that the previous one is OLS type
    if i == 1
        % Apply the SSC, filterbank and lastly the TNS
        AACSeq2(i,1).frameType = SSC(frameT,nextframeT,"OLS");
        AACSeq2(i,1).winType = window_type;
        frameF = filterbank(frameT,AACSeq2(i,1).frameType,AACSeq2(i,1).winType);
        if AACSeq2(i,1).frameType == "ESH"
            [AACSeq2(i,1).chl.frameF,AACSeq2(i,1).chl.TNScoeffs] = TNS(frameF(:,1:2:15),AACSeq2(i,1).frameType);
            [AACSeq2(i,1).chr.frameF,AACSeq2(i,1).chr.TNScoeffs] = TNS(frameF(:,2:2:16),AACSeq2(i,1).frameType);
        else
            [AACSeq2(i,1).chl.frameF,AACSeq2(i,1).chl.TNScoeffs] = TNS(frameF(:,1),AACSeq2(i,1).frameType);
            [AACSeq2(i,1).chr.frameF,AACSeq2(i,1).chr.TNScoeffs] = TNS(frameF(:,2),AACSeq2(i,1).frameType);
        end
        continue;
    end
    % Apply the SSC, filterbank and lastly the TNS
    AACSeq2(i,1).frameType = SSC(frameT,nextframeT,AACSeq2(i-1,1).frameType);
    AACSeq2(i,1).winType = window_type;

    % Apply the filterbank
    frameF = filterbank(frameT,AACSeq2(i,1).frameType,AACSeq2(i,1).winType);

    % If the frame type was ESH reshape the data to be 128-by-8 array
    % and assing to the corresponding channel
    if AACSeq2(i,1).frameType == "ESH"
        [AACSeq2(i,1).chl.frameF,AACSeq2(i,1).chl.TNScoeffs] = TNS(frameF(:,1:2:15),AACSeq2(i,1).frameType);
        [AACSeq2(i,1).chr.frameF,AACSeq2(i,1).chr.TNScoeffs] = TNS(frameF(:,2:2:16),AACSeq2(i,1).frameType);
    else
        [AACSeq2(i,1).chl.frameF,AACSeq2(i,1).chl.TNScoeffs] = TNS(frameF(:,1),AACSeq2(i,1).frameType);
        [AACSeq2(i,1).chr.frameF,AACSeq2(i,1).chr.TNScoeffs] = TNS(frameF(:,2),AACSeq2(i,1).frameType);
    end

end

end

