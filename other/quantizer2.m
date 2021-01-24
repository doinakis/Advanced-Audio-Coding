% a = -1:0.001:1;
% a = a';
% a = [a a];
% I = NaN(2001,2);
% minimum = -0.7;
% L = 2^4;
% delta = 0.1;
% I1  = ceil((a(a<=0)-minimum)./delta);
% I2  = floor((a(a>0)-minimum)./delta);
% I(a>0) = I2;
% I(a<=0) = I1;
% I(I >= L-1) = L-2;
% I(I<0) = 0;
% pq = minimum + I .* delta;
% plot(a,pq);
% hold on;
% plot(zeros(length(-1:0.1:1),1),-1:0.1:1,'k');
% plot(-1:0.1:1,zeros(length(-1:0.1:1),1),'k');
% title("4 bit symmetric quantizer with 0.1 step");
% yticks([-0.8:0.1:0.7])
% xticks([-0.85:0.1:0.75])

% 
% a = -1:0.001:1;
% a = a';
% minimum = -0.7;
% L = 2^4;
% delta = 0.1;
% I1  = ceil((a(a<=0)-minimum)./delta);
% I2  = floor((a(a>0)-minimum)./delta);
% I(a>0) = I2;
% I(a<=0) = I1;
% I(I >= L-1) = L-2;
% I(I<0) = 0;
% pq = minimum + I .* delta;
% plot(a,pq)
