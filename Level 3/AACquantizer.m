%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Function that implements the AAC quantizer. S matrix with the quantized
% symbols (1024-by-one), sfc matrix with the scalefactors (NB-by8 for ESH
% NB-by-1 otherwise), and G the gloabl gain of the current frame (1-by-8
% for ESH 1 otherwise).
% Where:
% frameF: The MDCT coefficients
% frameType: The type of the current frame
% SMR: The signal to Mask Ratio calculated for the corresponding frame
%%
function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)

global long_bands short_bands T
persistent MagicNumber MQ

if isempty(MQ) || isempty(MagicNumber)
    MagicNumber = 0.4054;
    MQ = 8191;
end

switch frameType
    case {"LSS","OLS","LPS"}
        % Calculate the acoustics threshold
        P = NaN(length(long_bands),1);
        for b = 1:length(long_bands)
            P(b,1) = sum(frameF(long_bands(b,2):long_bands(b,3)).^2);
        end
        T = P ./ SMR;

        % Initial estimation of the scalefactor gain
        a_hat = (16/3) * log2(max(frameF).^(3/4)/MQ);
        a = a_hat .* ones(length(long_bands),1);

        % Calculate the quantized symbol for every MDCT coefficient
        for b = 1:length(long_bands)
            a_all(long_bands(b,2):long_bands(b,3),1) = a(b);
        end
        S = sign(frameF) .* fix((abs(frameF) .* 2.^(-a_all/4)).^(3/4) + MagicNumber);

        % Check the power of quantization error
        X_hat = sign(S) .* (abs(S).^(4/3)) .* 2.^(a_all/4);
        Pe = NaN(length(long_bands),1);
        for b = 1:length(long_bands)
            Pe(b,1) = sum((frameF(long_bands(b,2):long_bands(b,3))- X_hat(long_bands(b,2):long_bands(b,3))).^2);
        end

        a(Pe <= T) = a(Pe <= T) + 1;
        % Variable thats indicates which a values should stop increasing
        stop_increasing = false(69,1);
        % If the quantization erro is less than the acoustics threshold
        while(any((Pe <= T) - stop_increasing))

            stop_increasing = stop_increasing | (Pe > T);

            for b = 1:length(long_bands)
                a_all(long_bands(b,2):long_bands(b,3),1) = a(b);
            end

            % Calculate the new symbols
            S = sign(frameF) .* fix((abs(frameF) .* 2.^(-a_all/4)).^(3/4) + MagicNumber);

            % Apply the inverse quantization
            X_hat = sign(S) .* (abs(S).^(4/3)) .* 2.^(a_all/4);

            % Calculate the power of the quantization error
            for b = 1:length(long_bands)
                Pe(b,1) = sum((frameF(long_bands(b,2):long_bands(b,3))-X_hat(long_bands(b,2):long_bands(b,3))).^2);
            end

            index = (Pe <= T) - stop_increasing;
            index(index<0) = 0;
            index = logical(index);
            a(index) = a(index) + 1;
            a(Pe > T) = a(Pe > T) - 1;
            % If the distance between a(b+1) - a(b) > 60 then stop increasing
            % the a values
            if max(abs(a(1:end-1) - a(2:end))) > 60
                a(index) = a(index) - 1;
                break;
            end
        end

        % Global gain
        G = fix(a(1));
        sfc = round(a(2:end) - a(1:end-1));

    otherwise
        P = NaN(length(short_bands),8);
        for b = 1:length(short_bands)
            P(b,:) = sum(frameF(short_bands(b,2):short_bands(b,3),:).^2,1);
        end

        T = P ./ SMR;

        % Initial estimation of the scalefactor gain
        a_hat = (16/3) * log2((max(frameF,[],1).^(3/4))/MQ);
        a = a_hat .* ones(length(short_bands),1);
        a_all = NaN(length(short_bands),8);
        for b = 1:length(short_bands)
            a_all(short_bands(b,2):short_bands(b,3),:) = repmat(a(b,:),length(short_bands(b,2):short_bands(b,3)),1);
        end
        S = sign(frameF) .* fix((abs(frameF) .* 2.^(-a_all/4)).^(3/4) + MagicNumber);

        % Check the power of quantization error
        X_hat = sign(S) .* (abs(S).^(4/3)) .* 2.^(a_all/4);
        Pe = NaN(length(short_bands),8);
        for b = 1:length(short_bands)
            Pe(b,:) = sum((frameF(short_bands(b,2):short_bands(b,3),:)- X_hat(short_bands(b,2):short_bands(b,3),:)).^2,1);
        end

        for i = 1:8

            a(Pe(:,i) <= T(:,i)) = a(Pe(:,i) <= T(:,i)) + 1;
            stop_increasing = false(42,1);

            while(any((Pe(:,i) <= T(:,i)) - stop_increasing))

                stop_increasing = stop_increasing | (Pe(:,i) > T(:,i));

                for b = 1:length(short_bands)
                    a_all(short_bands(b,2):short_bands(b,3),i) = repmat(a(b,i),length(short_bands(b,2):short_bands(b,3)),1);
                end

                % Calculate the quantized symbols
                S(:,i) = sign(frameF(:,i)) .* fix((abs(frameF(:,i)) .* 2.^(-a_all(:,i)/4)).^(3/4) + MagicNumber);

                % Apply the inverse quantization
                X_hat(:,i) = sign(S(:,i)) .* (abs(S(:,i)).^(4/3)) .* 2.^(a_all(:,i)/4);

                % Calculate the power of the quantization error
                for b = 1:length(short_bands)
                    Pe(b,i) = sum((frameF(short_bands(b,2):short_bands(b,3),i)- X_hat(short_bands(b,2):short_bands(b,3),i)).^2);
                end

                index = (Pe(:,i) <= T(:,i)) - stop_increasing;
                index(index<0) = 0;
                index = logical(index);
                a(index,i) = a(index,i) + 1;
                a(Pe(:,i) > T(:,i),i) = a(Pe(:,i) > T(:,i),i) - 1;
                % If the distance between a(b+1) - a(b) > 60 then we stop increasing
                % the a values
                if max(abs(a(1:end-1,i) - a(2:end,i))) > 60
                    a(index,i) = a(index,i) - 1;
                    break;
                end
            end
        end
        S = S(:);
        % Global gain
        G = fix(a(1,1:8));
        sfc = round(a(2:end,:) - a(1:end-1,:));

end

end

