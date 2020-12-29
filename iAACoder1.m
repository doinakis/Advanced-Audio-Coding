function x = iAACoder1(AACSeq1,fNameOut)

% Initialize empty signal x 
y = NaN(length(AACSeq1)*1024,2);

% Counter that helps reconstruct the initial signal
counter = 1025;

% For every frame in the struct
for i = 1:length(AACSeq1)
    
    % Initialize a frame with 2 channels
    frameF = [AACSeq1(i,1).chl.frameF AACSeq1(i,1).chr.frameF];
    
    if AACSeq1(i,1).frameType == "ESH"
        frameF(:,1:2:15) = AACSeq1(i,1).chl.frameF;
        frameF(:,2:2:16) = AACSeq1(i,1).chr.frameF;
    end
    % Apply the inverse filterbank to the corresponding frame
    frameT = ifilterbank(frameF,AACSeq1(i,1).frameType,AACSeq1(i,1).winType);
    
    % If its the first iteration then the singal x is set to be the whole
    % first frame
    if i == 1
        y(1:2048,:) = frameT;
        continue;
    end
    
    % If its not the first iteration then the first 1024 elements are added
    % to the previous frame's last 1024 elements  
    y(counter:counter+1023,:) = y(counter:counter+1023,:) + frameT(1:1024,:);
    counter = counter + 1024;
    
    % The other half of the frame is assinged to the next 1024 of the
    % singal
    y(counter:counter+1023,:) = frameT(1025:end,:);
end

% Write the audio file
audiowrite(fNameOut,y,48000);
% Return the signal values to variable x
if nargout == 1
    x = y;
end

global frameSNR frame1 frame2
frameSNR = zeros(length(AACSeq1),1);

for i = 1:length(AACSeq1)
    frameSNR(i) = snr([frame1(:,i) frame2(:,i)],x((1024*(i-1)+1):(1024*(i+1)),:) - [frame1(:,i) frame2(:,i)]); 
end
end

