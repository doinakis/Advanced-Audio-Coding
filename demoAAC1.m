function SNR = demoAAC1(fNameIn, fNameOut)

AACSeq1 = AACoder1(fNameIn);

y = audioread(fNameIn,'double');
x = iAACoder1(AACSeq1,fNameOut);


end

