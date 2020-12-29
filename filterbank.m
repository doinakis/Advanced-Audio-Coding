function frameF = filterbank(frameT,frameType,winType)

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

% Initialize an array for the mdct coefficients
frameF = NaN(N_long/2,2);


if (frameType == "ESH")
    % If the frame is ESH then calculate the 8 areas
    y = [];
    for i = 449:128:1345
        y = [y;frameT(i:i+255,:)];
    end
    % Apply the window to each area
    z = y .* window;
    % Apply the mdct
    counter_bottom = 1;
    counter_top = counter_bottom + 255;
    counter = 1;
    n = 0:N_short-1;
    n = n';
    n_zero = (N_short/2 + 1)/2;
    for i = 1:8
        for k = 0:N_short/2-1
            frameF(counter,:) = 2 .* sum(z(counter_bottom:counter_top,:) .* cos(2*pi/N_short * (n + n_zero)*(k + 1/2)));
            counter = counter + 1;
        end
        counter_bottom = counter_bottom + 256;
        counter_top = counter_top + 256;
    end
    
    temp1 = reshape(frameF(:,1),[128,8]);
    temp2 = reshape(frameF(:,2),[128,8]);
    frameF = [];
    for i = 1:8
        frameF = [frameF temp1(:,i) temp2(:,i)];
    end
else
    % Apply the window to the current frame
    z = frameT .* window;
    % Calculate the mdct 
    n = 0:N_long-1;
    n = n';
    n_zero = (N_long/2 + 1)/2;
    for k = 0:N_long/2-1
        frameF(k+1,:) = 2 .* sum(z .* cos(2*pi/N_long * (n + n_zero) * (k + 1/2)));
    end
end
end

