clc 
clear 
close all 
clear psycho
clear global
clear AACquantizer
clear iAACquantizer

tic
SNR = demoAAC3('../files/LicorDeCalandraca.wav', 'test','AACSeq');
toc

% x = NaN(69,69);

% global bval 
% bands = load('TableB219.mat');
% long_bands = bands.B219a;
% bval = long_bands(:,5);
% tic 
% for i = 1:69
%     for ii = 1:69
%         x(ii,i) = spreading_function(ii,i);
%     end
% end
% toc