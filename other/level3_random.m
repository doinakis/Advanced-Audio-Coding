clc
clear
close all


N = 2048;
% Hann window
n = (1:N)';
s = ones(N,1);
% 2 
sw = s .* (0.5 - 0.5 * cos((pi*(n+0.5))/(N/2)));

plot(1:2048,sw);

fast_fourier = fft(sw);
r = real(fast_fourier);
f = angle(fast_fourier);

r_pred = 2 * r_1 - r_2;
f_pred = 2 * f_1 - f_2;

