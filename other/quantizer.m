a = -1:0.001:1;
minimum = -0.7;
L = 2^4 - 1 ;
delta = 0.1;
I  = round((a-minimum)./delta);
I(I >= L) = L-1;
I(I<0) = 0;
pq = minimum + I .* delta;
plot(a,pq);
hold on;
plot(zeros(length(-1:0.1:1),1),-1:0.1:1,'k');
plot(-1:0.1:1,zeros(length(-1:0.1:1),1),'k');
yticks([-0.8:0.1:0.7])
xticks([-0.85:0.1:0.75])
xlabel("Value");
ylabel("Quantized value");
title("Quantizer");