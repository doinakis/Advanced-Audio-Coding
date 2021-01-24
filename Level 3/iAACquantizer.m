%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Function that implements the inverse AAC quantizer. And calculates the
% frameF (128-by-8 for ESH 1024-by-1 otherwise)
% Where:
% S: The quantized symbols
% sfc: The scalefactors of the bands
% frameType: The type of the corresponding frame
%%
function frameF = iAACquantizer(S, sfc, G, frameType)

global long_bands short_bands 

switch frameType
    case "ESH"
        
        % Reconstruct the a values
        a(1,:) = G + 10;
        for i = 2:length(sfc)
            a(i,:) = a(i-1,:) + sfc(i,:);
        end
        a_all = NaN(length(short_bands),8);
        for b = 1:length(short_bands)
            a_all(short_bands(b,2):short_bands(b,3),:) = repmat(a(b,:),length(short_bands(b,2):short_bands(b,3)),1);
        end
        S = reshape(S,128,8);
        % Reconstruct the frameF
        frameF = sign(S) .* (abs(S).^(4/3)) .* 2.^(a_all/4);
        
    otherwise
        
        a = NaN(length(sfc),1);
        % Reconstruct the a values
        a(1) = G + 10;
        for i = 2:length(sfc)
            a(i) = a(i-1) + sfc(i);
        end
        a_all = NaN(length(long_bands),1);
        for b = 1:length(long_bands)
            a_all(long_bands(b,2):long_bands(b,3),1) = a(b);
        end
        
        % Reconstruct the frameF
        frameF = sign(S) .* (abs(S).^(4/3)) .* 2.^(a_all/4);
end
end

