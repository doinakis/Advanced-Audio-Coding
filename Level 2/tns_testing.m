clc
clear
close all

 AACSeq1 = AACoder1("../files/LicorDeCalandraca.wav")

X = [AACSeq1(4).chl.frameF AACSeq1(4).chr.frameF];
counter = 1;
P = NaN(length(B219a),2);
for i = 1:length(B219a)
    if i == 69
        P(counter,:) = sum(X((B219a(i,2)):B219a(i,3),:).^2);
        continue;
    end
    P(counter,:) = sum(X((B219a(i,2)):B219a(i+1,2),:).^2);
    counter = counter + 1;
end

S = NaN(1024,2);
k = 1:1024;
for i = 1:69
    if i == 69
        lower = k >= B219a(i,2);
        upper = k <= B219a(i,3);
        index = find(lower == upper);
        S(index,:) = repmat(sqrt(P(i,:)),length(index),1);
        continue;
    end
    lower = k >= B219a(i,2);
    upper = k < B219a(i+1,2);
    index = find(lower == upper);
    S(index,:) = repmat(sqrt(P(i,:)),length(index),1);
end

for k = 1023:-1:1
    S(k,:) = (S(k,:) + S(k+1,:))./2;
    S(1025-k,:) = (S(1025-k,:) + S(1025-k-1,:))./2;
end
