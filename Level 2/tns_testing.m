clc
clear
close all

AACSeq1 = AACoder1("../files/LicorDeCalandraca.wav");
load("../files/TableB219.mat");
B219b(:,1:3) = B219b(:,1:3) + 1;

X = [AACSeq1(2).chl.frameF];
counter = 1;
P = NaN(length(B219b),8);
for i = 1:length(B219b)
    if i == length(B219b)
        P(counter,:) = sum(X((B219b(i,2)):B219b(i,3),:).^2);
        continue;
    end
    P(counter,:) = sum(X((B219b(i,2)):B219b(i+1,2),:).^2);
    counter = counter + 1;
end
    

S = NaN(128,8);
k = 1:128;
for i = 1:length(B219b)
    if i == length(B219b)
        lower = k >= B219b(i,2);
        upper = k <= B219b(i,3);
        index = find(lower == upper);
        S(index,:) = repmat(sqrt(P(i,:)),length(index),1);
        continue;
    end
    lower = k >= B219b(i,2);
    upper = k < B219b(i+1,2);
    index = find(lower == upper);
    S(index,:) = repmat(sqrt(P(i,:)),length(index),1);
end

for k = 127:-1:1
    S(k,:) = (S(k,:) + S(k+1,:))./2;
    S(129-k,:) = (S(129-k,:) + S(129-k-1,:))./2;
end

Xw = X./S;

a = NaN(4,8);
for i = 1:8
    [corellation1,lags] = xcorr(Xw(:,i),Xw(:,i),4,'biased');
    r = [corellation1(lags == 1);corellation1(lags == 2);
            corellation1(lags == 3);corellation1(lags == 4)];
    R = [corellation1(lags == 0) corellation1(lags == 1) corellation1(lags == 2) corellation1(lags == 3);
        corellation1(lags == -1) corellation1(lags == 0) corellation1(lags == 1) corellation1(lags == 2);
        corellation1(lags == -2) corellation1(lags == -1) corellation1(lags == 0) corellation1(lags == 1);
        corellation1(lags == -3) corellation1(lags == -2) corellation1(lags == -1) corellation1(lags == 0)];
    a(:,i) = R\r;
end

minimum = -0.85;
L = 2^4;
delta = 0.1;
I  = floor((a-minimum)./delta);
I(I >= L) = L-1;
I(I<0) = 0;
TNScoeffs = minimum + I .* delta;

for i = 1:8
    frameFout(:,i) = filter(TNScoeffs(:,i),1,X(:,i),[],1);
end

for i = 1:8
    ending(:,i) = filter(1,TNScoeffs(:,i),frameFout(:,i),[],1);
end

