function x = iAACoder1(AACSeq1,fNameOut)

frame1 = NaN(2048,length(AACSeq1));
frame2 = NaN(2048,length(AACSeq1));
x1 = NaN(1024,1);
x2 = NaN(1024,1);
for i = 1:length(AACSeq1)
    if AACSeq1(i).frameType == "ESH"
        frameF = [AACSeq1(i).chl.frameF AACSeq1(i).chr.frameF];
        frameF = reshape(frameF,[],1);
        frameT = ifilterbank(frameF,AACSeq1(i).frameType,AACSeq1(i).winType);
        frame1(:,i) = frameT(:,1);
        frame2(:,i) = frameT(:,2);
        if i == 1
            x1 = [x1;frame1(:,1)];
            x2 = [x2;frame2(:,2)];
            x1(i*1024:i*1024+1024) = x1(i*1024:i*1024+1024) + frame1(1:1024,i);
            x2(i*1024:i*1024+1024) = x2(i*1024:i*1024+1024) + frame2(1:1024,i);
            continue;
        end
        x1((((i-1)*1024)+1):((i-1)*1024)+1+1023) = x1((((i-1)*1024)+1):((i-1)*1024)+1+1023) + frame1(1:1024,i);
        x2((((i-1)*1024)+1):((i-1)*1024)+1+1023) = x2((((i-1)*1024)+1):((i-1)*1024)+1+1023) + frame2(1:1024,i);
        x1 = [x1;frame1(1025:end,1)];
        x2 = [x2;frame2(1025:end,2)];
    else
        frameT = ifilterbank([AACSeq1(i).chl.frameF AACSeq1(i).chr.frameF],AACSeq1(i).frameType,AACSeq1(i).winType);
        frame1(:,i) = frameT(:,1);
        frame2(:,i) = frameT(:,2);
        if i == 1
            x1 = [x1;frame1(:,1)];
            x2 = [x2;frame2(:,2)];
            continue;
        end
        x1((((i-1)*1024)+1):((i-1)*1024)+1+1023) = x1((((i-1)*1024)+1):((i-1)*1024)+1+1023) + frame1(1:1024,i);
        x2((((i-1)*1024)+1):((i-1)*1024)+1+1023) = x2((((i-1)*1024)+1):((i-1)*1024)+1+1023) + frame2(1:1024,i);
        x1 = [x1;frame1(1025:end,1)];
        x2 = [x2;frame2(1025:end,2)];
    end
end
x = [x1 x2];
audiowrite('attempt3.wav',x,48000);
end

