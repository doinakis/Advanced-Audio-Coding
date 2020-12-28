function x = iAACoder1(AACSeq1,fNameOut)

frame1 = NaN(2048,length(AACSeq1));
frame2 = NaN(2048,length(AACSeq1));
x = NaN(length(AACSeq1)*1024,2);
counter = 1025;
for i = 1:length(AACSeq1)
    frameF = [AACSeq1(i,1).chl.frameF AACSeq1(i,1).chr.frameF];
    if AACSeq1(i,1).frameType == "ESH"
        frameF = [reshape(AACSeq1(i,1).chl.frameF,[],1) reshape(AACSeq1(i,1).chl.frameF,[],1)];
    end
    frameT = ifilterbank(frameF,AACSeq1(i,1).frameType,AACSeq1(i).winType);
    
    if i == 1
        x(1:2048,:) = frameT;
        continue;
    end
    
    x(counter:counter+1023,:) = x(counter:counter+1023,:) + frameT(1:1024,:);
    counter = counter + 1024;
    x(counter:counter+1023,:) = frameT(1025:end,:);
end
audiowrite(fNameOut,x,48000);
end

