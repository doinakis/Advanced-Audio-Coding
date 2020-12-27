function AACSeq1 = AACoder1(fNameIn)

[y,Fs] = audioread(fNameIn,'double');

% Frame Creation 
frame_counter = 1;
for i = 1:1024:size(y,1)
    if i == 1 
        frame1(:,frame_counter)=[zeros(1024,1);y(1:1024,1)];
        frame2(:,frame_counter)=[zeros(1024,1);y(1:1024,2)];
        frame_counter = frame_counter + 1;
        continue;
    end
    if i+2047 > size(y,1)
        frame1(1:2048,frame_counter) = [y(i:end,1);zeros(i+2047-size(y,1),1)];
        frame2(1:2048,frame_counter) = [y(i:end,2);zeros(i+2047-size(y,1),1)];
        break;
    end
    frame1(1:2048,frame_counter) = y(i:i+2047,1);
    frame2(1:2048,frame_counter) = y(i:i+2047,2);
    frame_counter = frame_counter + 1;
end

for i = 1:276
    frameT = [frame1(:,i) frame2(:,i)];
    
    if i == 276 
        nextframeT = zeros(2048,2);
    else
        nextframeT = [frame1(:,i+1) frame2(:,i+1)];
    end
    if i == 1
        prevframeType = "OLS";
        AACSeq1(i).frameType = SSC(frameT,nextframeT,prevframeType);
        AACSeq1(i).winType = "SIN";
        frameF = filterbank(frameT,AACSeq1(i).frameType,AACSeq1(i).winType);
        if AACSeq1(i).frameType == "ESH"
            AACSeq1(i).chl.frameF = reshape(frameF(:,1),[128,8]);
            AACSeq1(i).chr.frameF = reshape(frameF(:,2),[128,8]);
        else
            AACSeq1(i).chl.frameF = frameF(:,1);
            AACSeq1(i).chr.frameF = frameF(:,2);
        end
        continue;
    end
    prevframeType = AACSeq1(i-1).frameType;
    AACSeq1(i).frameType = SSC(frameT,nextframeT,prevframeType);
    AACSeq1(i).winType = "SIN";
    frameF = filterbank(frameT,AACSeq1(i).frameType,AACSeq1(i).winType);
    
    if AACSeq1(i).frameType == "ESH"
        AACSeq1(i).chl.frameF = reshape(frameF(:,1),[128,8]);
        AACSeq1(i).chr.frameF = reshape(frameF(:,2),[128,8]);
    else
        AACSeq1(i).chl.frameF = frameF(:,1);
        AACSeq1(i).chr.frameF = frameF(:,2);
    end

end

end

