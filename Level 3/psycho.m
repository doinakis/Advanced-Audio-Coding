%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Function that calculates the Signal to Mask Ratio. If the frame is ESH
% the output is 42-by-8, otherwise it is 69-by-1.
% Where:
% frameT: The current frame
% frameType: The type of the current frame 
% frameTprev1: The values of the previous frame
% frameTprev2: The values of the frameTprev1's previous frame
%%
function SMR = psycho(frameT, frameType, frameTprev1, frameTprev2)

% Create persistent variables for efficiency
persistent spreading_func_ESH spreading_func_OTHER
persistent q_hat_long q_hat_short hann_long hann_short

global long_bands short_bands

if isempty(spreading_func_ESH) || isempty(spreading_func_OTHER)
    qsthr_long =  long_bands(:,6);
    qsthr_short = short_bands(:,6);
    q_hat_long = eps * 1024 * 10.^(qsthr_long ./ 10);
    q_hat_short = eps * 128 * 10.^(qsthr_short ./ 10);
    bval = long_bands(:,5);
    for ii = 1:length(bval)
        for i = 1:length(bval)
            spreading_func_OTHER(ii,i) = spreading_function(ii,i,bval);
        end
    end
    bval = short_bands(:,5);
    for ii = 1:length(bval)
        for i = 1:length(bval)
            spreading_func_ESH(ii,i) = spreading_function(ii,i,bval);
        end
    end
    N = 256;
    hann_short = 0.5 - 0.5 * cos(2*pi/N * ((0:N-1)' + 0.5));
    hann_short = repmat(hann_short,8,1);
    N = 2048;
    hann_long = 0.5 - 0.5 * cos(2*pi/N * ((0:N-1)' + 0.5));
end

switch frameType
    case "ESH"
        % N
        N = 256;
        
        % If the frame is ESH then calculate the 8 areas
        s = NaN(2048,1);
        prev_Frames = NaN(2048,1);
        
        % The y variable holds 2048-by-2 samples. It holds the 8*256 subframes
        counter = 1;
        for i = 449:128:1345
            s(counter:counter+255,:) = frameT(i:i+255,:);
            prev_Frames(counter:counter+255,:) = frameTprev1(i:i+255,:);
            counter = counter + 256;
        end
        
        % Apply the filter to every subframe 
        sw = s .* hann_short;
        prev_Frames = prev_Frames .* hann_short;
        sw = reshape(sw,[N,8]);
        prev_Frames = reshape(prev_Frames,[N,8]);
        
        % Only the 7th and 8th subframes are needed from the previous frame 
        prev_Frames = prev_Frames(:,7:8);
        
        % Calculate the fft of every sub-frame
        fft_sw = fft(sw,[],1);
        fft_sw = fft_sw(1:128,:);
        rw = abs(fft_sw);
        fw = angle(fft_sw);
        fft_prev = fft(prev_Frames);
        fft_prev = fft_prev(1:128,:);
        r_1 = abs(fft_prev(:,2));
        r_2 = abs(fft_prev(:,1));
        f_1 = angle(fft_prev(:,2));
        f_2 = angle(fft_prev(:,1));
        r_pred = NaN(128,8);
        f_pred = NaN(128,8);
        
        % Calculate the predicted values for every sub-frame
        for i = 1:8 
            r_pred(:,i) = 2 .* r_1 - r_2;
            f_pred(:,i) = 2 .* f_1 - f_2;
            r_2 = r_1;
            f_2 = f_1;
            r_1 = rw(:,i);
            f_1 = fw(:,i);
        end
        
        % Calculate the predictabilities
        cw = sqrt((rw.*cos(fw) - r_pred.*cos(f_pred)).^2 + (rw.*sin(fw) - r_pred.*sin(f_pred)).^2)./(rw + abs(r_pred));
        
        e = NaN(length(short_bands),8);
        c = NaN(length(short_bands),8);
        % 
        for b = 1:length(short_bands)
            e(b,:) = sum(rw(short_bands(b,2):short_bands(b,3),:).^2,1);
            c(b,:) = sum(cw(short_bands(b,2):short_bands(b,3),:).*rw(short_bands(b,2):short_bands(b,3),:).^2,1);
        end
        
        ecb = NaN(length(short_bands),8);
        ct = NaN(length(short_bands),8);
        % Combine the energy and predictability with the spreading function
        for b = 1:length(short_bands)
            ecb(b,:) = sum(e .* spreading_func_ESH(:,b),1);
            ct(b,:) = sum(c .* spreading_func_ESH(:,b),1);
        end
        
        % Normalization
        cb = ct ./ ecb;
        spreading_sum = sum(spreading_func_ESH,1)';
        en = ecb ./ spreading_sum;
        
        % Calculate tonality index
        tb = -0.299 - 0.43 .* log(cb);
        
        % tb values are in the interval (0,1)
        tb = (tb - min(tb))./(max(tb) - min(tb));
        
        % Noise Masking Tone = 6dB
        NMT = 6;
       
        % Tone Masking Noise = 18dB
        TMN = 18;
        
        %SNR calculation for every b 
        SNR = tb .* TMN + (1 - tb) .* NMT; 
        
        % Convert db to energy
        bc = 10.^(-SNR./10);
        
        % Calculate the energy threshold
        nb = en .* bc;
        
        % Calculate the noise level for every band
        npart = max(nb,q_hat_short);
        
        % Calculate Signal to Mask Ratio
        SMR = e ./ npart;
            
    otherwise
        
        % Apply the window to the current frame calculate the frame's fft
        % and keep only the first 1024 samples. From the fft extract the
        % absolute value and the phase for each sample
        sw = frameT .* hann_long;
        fft_sw = fft(sw);
        fft_sw = fft_sw(1:1024);
        rw = abs(fft_sw);
        fw = angle(fft_sw);
        sw_1 = frameTprev1 .* hann_long;
        sw_2 = frameTprev2 .* hann_long;
        fft_sw1 = fft(sw_1);
        fft_sw1 = fft_sw1(1:1024);
        rw_1 = abs(fft_sw1);
        fw_1 = angle(fft_sw1);
        fft_sw2 = fft(sw_2);
        fft_sw2 = fft_sw2(1:1024);
        rw_2 = abs(fft_sw2);
        fw_2 = angle(fft_sw2);
        
        
        % Calculate the predicted values 
        r_pred = 2 .* rw_1 - rw_2;
        f_pred = 2 .* fw_1 - fw_2;
        
        % Calculate the predictabilities 
        cw = sqrt((rw.*cos(fw) - r_pred.*cos(f_pred)).^2 + (rw.*sin(fw) - r_pred.*sin(f_pred)).^2)./(rw + abs(r_pred));
        
        e = NaN(length(long_bands),1);
        c = NaN(length(long_bands),1);
        
        %  
        for b = 1:length(long_bands)
            e(b,1) = sum(rw(long_bands(b,2):long_bands(b,3)).^2);
            c(b,1) = sum(cw(long_bands(b,2):long_bands(b,3)).*rw(long_bands(b,2):long_bands(b,3)).^2);
        end
        
        ecb = NaN(length(long_bands),1);
        ct = NaN(length(long_bands),1);
        
        % Combine the energy and predictability with the spreading function 
        for b = 1:length(long_bands)
            ecb(b,1) = sum(e .* spreading_func_OTHER(:,b));
            ct(b,1) = sum(c .* spreading_func_OTHER(:,b));
        end
        
        % Normalization 
        cb = ct ./ ecb;
        spreading_sum = sum(spreading_func_OTHER,1)';
        en = ecb ./ spreading_sum;
        
        % Calculate the tonality indices
        tb = -0.299 - 0.43 .* log(cb);
        
        % tb values are in the interval (0,1)
        tb = (tb - min(tb))./(max(tb) - min(tb));
        
        % Noise Masking Tone = 6dB
        NMT = 6;
        
        % Tone Masking Noise = 18dB
        TMN = 18;
        
        %SNR calculation for every b
        SNR = tb .* TMN + (1 - tb) .* NMT;
        
        % Convert db to energy
        bc = 10.^(-SNR./10);
        
        % Calculate the energy threshold
        nb = en .* bc;
        
        % Calculate the noise level for every band
        npart = max(nb,q_hat_long);
        
        % Calculate Signal to Mask Ratio
        SMR = e ./ npart;
       
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that implements the spreading function calculation
% Where:
% i: the band that causes the spreading 
% j: the band that the spreading affects
%%
function x = spreading_function(i,j,bval)
if i>=j
    tmpx = 3 * (bval(j)-bval(i));
else
    tmpx = 1.5 * (bval(j)-bval(i));
end

tmpz = 8 * min((tmpx-0.5)^2 - 2*(tmpx-0.5),0);
tmpy = 15.811389 + 7.5 * (tmpx+0.474) - 17.5 * (1+(tmpx+0.474)^2)^(1/2);

if tmpy<-100
    x = 0;
else
    x = 10^((tmpz+tmpy)/10);
end

end

