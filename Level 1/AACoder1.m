%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that calculates the coded frames and returns a K-by-1 srtuct
% with the frameType,winType,chl.frameF (the left channels mdct coefficients),
% chr.frameF (the right channels mdct coefficients). If the frame is ESH
% then the chr and chl contains a 128-by-8 matrix, otherwise contain a
% 1024-by-1. Where:
% fNameIn: The path (or the name if its in the same folder) of the .wav
% file
%%
function AACSeq1 = AACoder1(fNameIn)

% Read the audio signal from the input
y = audioread(fNameIn,'double');

% The window type 
window_type = "SIN";

% Preallocate the space for the frames 
AACSeq1(ceil(size(y,1)/1024),1) = struct();

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
        AACSeq1(i,1).frameType = SSC(frameT,nextframeT,"OLS");
        AACSeq1(i,1).winType = window_type;
        frameF = filterbank(frameT,AACSeq1(i,1).frameType,AACSeq1(i,1).winType);
        if AACSeq1(i,1).frameType == "ESH"
            AACSeq1(i,1).chl.frameF = frameF(:,1:2:15);
            AACSeq1(i,1).chr.frameF = frameF(:,2:2:16);
        else
            AACSeq1(i,1).chl.frameF = frameF(:,1);
            AACSeq1(i,1).chr.frameF = frameF(:,2);
        end
        continue;
    end
    
    AACSeq1(i,1).frameType = SSC(frameT,nextframeT,AACSeq1(i-1,1).frameType);
    AACSeq1(i,1).winType = window_type;
    
    % Apply the filterbank
    frameF = filterbank(frameT,AACSeq1(i,1).frameType,AACSeq1(i,1).winType);
    
    % If the frame type was ESH reshape the data to be 128-by-8 array
    % and assing to the corresponding channel
    if AACSeq1(i,1).frameType == "ESH"
        AACSeq1(i,1).chl.frameF = frameF(:,1:2:15);
        AACSeq1(i,1).chr.frameF = frameF(:,2:2:16);
    else
        AACSeq1(i,1).chl.frameF = frameF(:,1);
        AACSeq1(i,1).chr.frameF = frameF(:,2);
    end

end

end

