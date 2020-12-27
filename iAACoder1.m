function x = iAACoder1(AACSeq1,fNameOut)

frame1 = NaN(2048,length(AACSeq1));
frame2 = NaN(2048,length(AACSeq1));
x1 = [];
x2 = [];
for i = 1:length(AACSeq1)
    if AACSeq1(i).frameType == "ESH"
        frameF = [AACSeq1(i).chl.frameF AACSeq1(i).chr.frameF];
        frameF = reshape(frameF,[],1);
        frameT = ifilterbank(frameF,AACSeq1(i).frameType,AACSeq1(i).winType);
        frame1(:,i) = frameT(:,1);
        frame2(:,i) = frameT(:,2);
        x1 = [x1;frameT(1:1024,1)];
        x2 = [x2;frameT(1:1024,2)];
        
    else
        frameT = ifilterbank([AACSeq1(i).chl.frameF AACSeq1(i).chr.frameF],AACSeq1(i).frameType,AACSeq1(i).winType);
        frame1(:,i) = frameT(:,1);
        frame2(:,i) = frameT(:,2);
        x1 = [x1;frameT(1:1024,1)];
        x2 = [x2;frameT(1:1024,2)];
        
    end
end
x = [x1 x2];
end

