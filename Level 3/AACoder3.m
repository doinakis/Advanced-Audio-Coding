%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------- Implements Level 3 of the Assignment (PSYCHO) --------------
% Function that calculates the coded frames for the compressed signal
% and returns a K-by-1 srtuct with the frameType,winType,chl.TNScoeffs,chr.TNScoeffs
% (the quantized TNS coefficients of the left and right channel, 4-by-8 for ESH 4-by-1
% else), chl.T,chr.T (the left and right channel's acoustics threshold),
% chl.G,chr.G (the left and right channel's global gain, 8 for ESH 1
% else),chl.sfc,chr.sfc (the left and right channel's encoded sequence
% sfc),chl.stream,chr.stream (the left and right channel's encoded
% quantized MDCT sequence), chl.codebook, chr.codebook (the left and right
% channel's codebook)
% Where:
% fNameIn: The path (or the name if its in the same folder) of the .wav
% file
% fnameAACoded: The desirable name of the matrix that will hold the
% compressed sigal
%%
function AACSeq3 = AACoder3(fNameIn, fnameAACoded)

global long_bands short_bands T

% Load the Psychoacoustic Model bands
bands = load('TableB219.mat');
long_bands = bands.B219a;
short_bands = bands.B219b;
long_bands(:,1:3) = long_bands(:,1:3) + 1;
short_bands(:,1:3) = short_bands(:,1:3) + 1;

huffLUT = loadLUT();
% Read the audio signal from the input
y = audioread(fNameIn,'double');

% The window type
window_type = "SIN";

% Preallocate the space for the frames
AACSeq3(ceil(size(y,1)/1024),1) = struct();

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
% Initialize the previous frames
prevframe1 = zeros(2048,2);
prevframe2 = zeros(2048,2);

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

    % Create dual channel previous frames
    % For the first frame assume that the previous one is OLS type
    if i == 1
        % Apply the SSC, filterbank and lastly the TNS
        AACSeq3(i,1).frameType = SSC(frameT,nextframeT,"OLS");
        AACSeq3(i,1).winType = window_type;
        frameF = filterbank(frameT,AACSeq3(i,1).frameType,AACSeq3(i,1).winType);
        if AACSeq3(i,1).frameType == "ESH"
            [frameF_left,AACSeq3(i,1).chl.TNScoeffs] = TNS(frameF(:,1:2:15),AACSeq3(i,1).frameType);
            [frameF_right,AACSeq3(i,1).chr.TNScoeffs] = TNS(frameF(:,2:2:16),AACSeq3(i,1).frameType);
        else
            [frameF_left,AACSeq3(i,1).chl.TNScoeffs] = TNS(frameF(:,1),AACSeq3(i,1).frameType);
            [frameF_right,AACSeq3(i,1).chr.TNScoeffs] = TNS(frameF(:,2),AACSeq3(i,1).frameType);
        end
        % Apply the psychoacoustic model for both channels 
        SMR_left = psycho(frameT(:,1), AACSeq3(i,1).frameType, prevframe1(:,1), prevframe2(:,1));
        SMR_right = psycho(frameT(:,2), AACSeq3(i,1).frameType, prevframe1(:,2), prevframe2(:,2));
        
        % Apply the quantization based on the acoustics threshold
        [S_left, AACSeq3(i,1).chl.sfc, AACSeq3(i,1).chl.G] = AACquantizer(frameF_left, AACSeq3(i,1).frameType, SMR_left);
        AACSeq3(i,1).chl.T = T;
        [S_right, AACSeq3(i,1).chr.sfc, AACSeq3(i,1).chr.G] = AACquantizer(frameF_right, AACSeq3(i,1).frameType, SMR_right);
        AACSeq3(i,1).chr.T = T;
        
        % Encode both streams and sfc coefficients with huffman
        [AACSeq3(i,1).chl.stream, AACSeq3(i,1).chl.codebook] = encodeHuff(S_left, huffLUT);
        [AACSeq3(i,1).chr.stream, AACSeq3(i,1).chr.codebook] = encodeHuff(S_right, huffLUT);
        AACSeq3(i,1).chl.sfc = encodeHuff(AACSeq3(i,1).chl.sfc(:), huffLUT, 12);
        AACSeq3(i,1).chr.sfc = encodeHuff(AACSeq3(i,1).chr.sfc(:), huffLUT, 12);
        
        % Set the prev frames for the next iteration
        prevframe2 = prevframe1;
        prevframe1 = [frame1(:,i) frame2(:,i)];
        continue;
    end
    % Apply the SSC, filterbank and lastly the TNS
    AACSeq3(i,1).frameType = SSC(frameT,nextframeT,AACSeq3(i-1,1).frameType);
    AACSeq3(i,1).winType = window_type;

    % Apply the filterbank
    frameF = filterbank(frameT,AACSeq3(i,1).frameType,AACSeq3(i,1).winType);

    % If the frame type was ESH reshape the data to be 128-by-8 array
    % and assing to the corresponding channel
    if AACSeq3(i,1).frameType == "ESH"
        [frameF_left,AACSeq3(i,1).chl.TNScoeffs] = TNS(frameF(:,1:2:15),AACSeq3(i,1).frameType);
        [frameF_right,AACSeq3(i,1).chr.TNScoeffs] = TNS(frameF(:,2:2:16),AACSeq3(i,1).frameType);
    else
        [frameF_left,AACSeq3(i,1).chl.TNScoeffs] = TNS(frameF(:,1),AACSeq3(i,1).frameType);
        [frameF_right,AACSeq3(i,1).chr.TNScoeffs] = TNS(frameF(:,2),AACSeq3(i,1).frameType);
    end
    
    % Apply the psychoacoustic model for both channels 
    SMR_left = psycho(frameT(:,1), AACSeq3(i,1).frameType, prevframe1(:,1), prevframe2(:,1));
    SMR_right = psycho(frameT(:,2), AACSeq3(i,1).frameType, prevframe1(:,2), prevframe2(:,2));
    
    % Apply the quantization based on the acoustics threshold
    [S_left, AACSeq3(i,1).chl.sfc, AACSeq3(i,1).chl.G] = AACquantizer(frameF_left, AACSeq3(i,1).frameType, SMR_left);
    AACSeq3(i,1).chl.T = T;
    [S_right, AACSeq3(i,1).chr.sfc, AACSeq3(i,1).chr.G] = AACquantizer(frameF_right, AACSeq3(i,1).frameType, SMR_right);
    AACSeq3(i,1).chr.T = T;
    
    % Encode both streams and sfc coefficients with huffman
    [AACSeq3(i,1).chl.stream, AACSeq3(i,1).chl.codebook] = encodeHuff(S_left, huffLUT);
    [AACSeq3(i,1).chr.stream, AACSeq3(i,1).chr.codebook] = encodeHuff(S_right, huffLUT);
    AACSeq3(i,1).chl.sfc = encodeHuff(AACSeq3(i,1).chl.sfc(:), huffLUT, 12);
    AACSeq3(i,1).chr.sfc = encodeHuff(AACSeq3(i,1).chr.sfc(:), huffLUT, 12);

    % Set the prev frames for the next iteration
    prevframe2 = prevframe1;
    prevframe1 = [frame1(:,i) frame2(:,i)];
end

save(fnameAACoded,'AACSeq3');
end

