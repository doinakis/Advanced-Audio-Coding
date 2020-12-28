function frameT = ifilterbank(frameF,frameType,winType)

% Define the sizes of long and short windows
N_long = 2048;
N_short = 256;

switch winType
    % Create with Kaiser Bessel Derived
    case "KBD"
        % Create Wl window
        win_left = NaN(N_long/2,1);
        win_right = NaN(N_long/2,1);
        a = 4;
        % Compute the w(n) and calculate the sum of its values
        w = kaiser(N_long/2 + 1,pi*a);
        w_sum = sum(w);
        % Create the left and the rigth window
        for n = 1:N_long/2
            win_left(n) = sqrt(sum(w(1:n))/w_sum);
            win_right(n) = sqrt(sum(w(N_long/2+1-n:-1:1))/w_sum);
        end
        win_long = [win_left;win_right];
        % Create Ws window
        win_left = NaN(N_short/2,1);
        win_right = NaN(N_short/2,1);
        a = 6;
        w = kaiser(N_short/2+1,pi*a);
        w_sum = sum(w);
        for n = 1:N_short/2
            win_left(n) = sqrt(sum(w(1:n))/w_sum);
            win_right(n) = sqrt(sum(w(N_short/2+1-n:-1:1))/w_sum);
        end
        win_short = [win_left;win_right];
    % Create the sin windows
    case "SIN"
        win_short = sin(pi/N_short * ((0:N_short-1) + 1/2))';
        win_long = sin(pi/N_long * ((0:N_long-1) + 1/2))';
end

% Create the window according to the frame type
switch frameType
    case "OLS"
        window = win_long;
    case "ESH"
        window = win_short;
        window = repmat(window,8,1);
    case "LSS"
        window = [win_long(1:1024);ones(448,1);win_short(129:end);zeros(448,1)];
    case "LPS"
        window = [zeros(448,1);win_short(1:128);ones(448,1);win_long(1025:end)];
end

frameT = NaN(N_long,2);

if (frameType == "ESH")
    counter_bottom = 1;
    counter_top = counter_bottom + 127;
    counter = 1;
    k = 0:N_short/2-1;
    k = k';
    n_zero = (N_short/2 + 1)/2;
    for i = 1:8
        for n = 0:N_short-1
            frameT(counter,:) = 2/N_short .* sum(frameF(counter_bottom:counter_top,:) .* cos(2*pi/N_short * (n + n_zero)*(k + 1/2)));
            counter = counter + 1;
        end
        counter_bottom = counter_bottom + 128;
        counter_top = counter_top + 128;
    end
    frameT = frameT .*window;
    channel1 = reshape(frameT(:,1),[256,8]);
    channel2 = reshape(frameT(:,2),[256,8]);
    total1 = zeros(448,1);
    total2 = zeros(448,1);
    total1 = [total1;channel1(:,1)];
    total2 = [total2;channel2(:,2)];

    for i = 2:8
        total1 = [total1; channel1(129:129 + 127,i-1) + channel1(1:128,i)];
        total2 = [total2; channel2(129:129 + 127,i-1) + channel2(1:128,i)];
    end
    
    total1 = [total1;zeros(448,1)];
    total2 = [total2;zeros(448,1)];
    frameT = [total1 total2];
    
else

    n_zero = (N_long/2 + 1)/2;
    k = 0:N_long/2-1;
    k = k';
    for n = 0:N_long-1
        frameT(n+1,:) = 2/N_long .* sum(frameF .* cos(2*pi/N_long * (n + n_zero)*(k + 1/2)));
    end
    frameT = frameT .* window;
end

end