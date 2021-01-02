function [frameFout,TNScoeffs] = TNS(frameFin,frameType)

bands = load('TableB219.mat');
long_bands = bands.B219a;
short_bands = bands.B219b;
long_bands(:,1:3) = long_bands(:,1:3) + 1;
short_bands(:,1:3) = short_bands(:,1:3) + 1;

switch frameType 
    case "ESH"
        bands = short_bands;
        n = 128;
        z = 8;
        frameFout = NaN(n,8);
    otherwise
        bands = long_bands;
        z = 1;
        n = 1024;
        frameFout = NaN(n,1);
end



counter = 1;
P = NaN(length(bands),z);
for i = 1:length(bands)
    if i == length(bands)
        P(counter,:) = sum(frameFin((bands(i,2)):bands(i,3),:).^2);
        continue;
    end
    P(counter,:) = sum(frameFin((bands(i,2)):bands(i+1,2),:).^2);
    counter = counter + 1;
end

S = NaN(n,z);
k = 1:n;
for i = 1:length(bands)
    if i == length(bands)
        lower = k >= bands(i,2);
        upper = k <= bands(i,3);
        index = find(lower == upper);
        S(index,:) = repmat(sqrt(P(i,:)),length(index),1);
        continue;
    end
    lower = k >= bands(i,2);
    upper = k < bands(i+1,2);
    index = find(lower == upper);
    S(index,:) = repmat(sqrt(P(i,:)),length(index),1);
end

for k = n-1:-1:1
    S(k,:) = (S(k,:) + S(k+1,:))./2;
    S(n+1-k,:) = (S(n+1-k,:) + S(n+1-k-1,:))./2;
end

Xw = frameFin./S;

if frameType == "ESH"
    a = NaN(4,8);
    for i = 1:8
        [corellation,lag] = xcorr(Xw(:,i),Xw(:,i),4,'biased');
        r = [corellation(lag == 1);corellation(lag == 2);
            corellation(lag == 3);corellation(lag == 4)];
        R = [corellation(lag == 0) corellation(lag == 1) corellation(lag == 2) corellation(lag == 3);
            corellation(lag == -1) corellation(lag == 0) corellation(lag == 1) corellation(lag == 2);
            corellation(lag == -2) corellation(lag == -1) corellation(lag == 0) corellation(lag == 1);
            corellation(lag == -3) corellation(lag == -2) corellation(lag == -1) corellation(lag == 0)];
        a(:,i) = R\r;
    end
else
    
    [corellation,lag] = xcorr(Xw,Xw,4,'biased');
    r = [corellation(lag == 1);corellation(lag == 2);
        corellation(lag == 3);corellation(lag == 4)];
    R = [corellation(lag == 0) corellation(lag == 1) corellation(lag == 2) corellation(lag == 3);
        corellation(lag == -1) corellation(lag == 0) corellation(lag == 1) corellation(lag == 2);
        corellation(lag == -2) corellation(lag == -1) corellation(lag == 0) corellation(lag == 1);
        corellation(lag == -3) corellation(lag == -2) corellation(lag == -1) corellation(lag == 0)];
    a = R\r;
end

% Quantization of a 
minimum = -0.85;
L = 2^4;
delta = 0.1;
I  = floor((a-minimum)./delta);
I(I >= L) = L-1;
I(I<0) = 0;
TNScoeffs = minimum + I .* delta;
TNScoeffs = -TNScoeffs;

if frameType == "ESH"
    for i = 1:8
        frameFout(:,i) = filter([1;TNScoeffs(:,i)],1,frameFin(:,i),[],1);
    end
else
    frameFout = filter([1;TNScoeffs],1,frameFin,[],1);
end

end

