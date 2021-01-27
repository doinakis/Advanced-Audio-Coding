clc 
clear 
close all 
clear psycho
clear global
clear AACquantizer
clear iAACquantizer

tic
[SNR, bitrate, compression] = demoAAC3('../files/LicorDeCalandraca.wav', 'test','AACSeq');
toc