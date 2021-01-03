%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that calculates the MDCT coeffs after TNS is applied.(frameFout)
% Where:
% frameFin: The MDCT coefficients 
% frameType: The type of the input frame
%%
function [frameFout,TNScoeffs] = TNS(frameFin,frameType)

% Load the Psychoacoustic Model bands
bands = load('TableB219.mat');
long_bands = bands.B219a;
short_bands = bands.B219b;
long_bands(:,1:3) = long_bands(:,1:3) + 1;
short_bands(:,1:3) = short_bands(:,1:3) + 1;

switch frameType 
    case "ESH"
        % If the frame is ESH then use short bands
        bands = short_bands;
        n = 128;
        z = 8;
        frameFout = NaN(n,8);
    otherwise
        % For every other type use the long bands
        bands = long_bands;
        z = 1;
        n = 1024;
        frameFout = NaN(n,1);
end

% Calculate the energy of every band
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

% Calculate the normalization coefficients S
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

% Smoothing of the coefficients 
for k = n-1:-1:1
    S(k,:) = (S(k,:) + S(k+1,:))./2;
    S(n+1-k,:) = (S(n+1-k,:) + S(n+1-k-1,:))./2;
end

% Normalization of the MDCT coefficients acording to the band they belong
% to
Xw = frameFin./S;

% Calculate the coefficients that miniizes the mean square error
if frameType == "ESH"
    % If the band is ESH then it calculates the coefficients for every
    % sub-frame
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

% Quantization of the coefficients so that the inverse filter can be
% applied 
% Quantization using R = 4 bits quantizer with 0.1 step
minimum = -0.8;
L = 2^4;
delta = 0.1;
I  = round((a-minimum)./delta);
I(I >= L) = L-1;
I(I<0) = 0;
TNScoeffs = minimum + I .* delta;

% Apply the filter H(z) = 1-a(1)*z^(-1) - ... - a(4)*z^(-4).
if frameType == "ESH"
    for i = 1:8
        if ~isstable(1,[1;-TNScoeffs(:,i)])
            disp('ERROR: Filter not reversable');
            exit(1);
        end
        frameFout(:,i) = filter([1;-TNScoeffs(:,i)],1,frameFin(:,i),[],1);
    end
else
    if ~isstable(1,[1;-TNScoeffs])
        disp('ERROR: Filter not reversable');
        exit(1);
    end
    frameFout = filter([1;-TNScoeffs],1,frameFin,[],1);
end

end

