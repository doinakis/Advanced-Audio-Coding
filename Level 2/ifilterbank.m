%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that implements the inverse filterbank. It returns the frame
% in time with coefficients frameF. It contains a 2048-by-2 array. If the
% frame its ESH the first and last 448 elements ara zero padded.
% Where:
% frameT: The frame in time
% frameType: The type of the frame (OLS,LPS,LSS,ESH)
% winType: The type of the window to be applied (SIN,KBD)
%%
function frameT = ifilterbank(frameF,frameType,winType)

% MDCT matrices coefficients
persistent W_long W_short

if isempty(W_long)
    % Calculate the iMDCT coefficients matrices for short and long windows
    N_long = 2048;
    n_zero = (N_long/2 + 1)/2;
    k = 0:N_long/2-1;
    n = (0:N_long-1)';
    W_long = 2/N_long * cos(2*pi/N_long .* (n + n_zero) .* (k + 1/2));
    
    N_short = 256;
    n_zero = (N_short/2 + 1)/2;
    k = 0:N_short/2-1;
    n = (0:N_short-1)';
    W_short = 2/N_short * cos(2*pi/N_short .* (n + n_zero) .* (k + 1/2));
end

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

        % Combine the left and right windows in one long window
        win_long = [win_left;win_right];

        % Create Ws window
        win_left = NaN(N_short/2,1);
        win_right = NaN(N_short/2,1);
        a = 6;
        w = kaiser(N_short/2 +1 ,pi*a);
        w_sum = sum(w);

        for n = 1:N_short/2
            win_left(n) = sqrt(sum(w(1:n))/w_sum);
            win_right(n) = sqrt(sum(w(N_short/2+1-n:-1:1))/w_sum);
        end

        % Combine the left and right windows in one short window
        win_short = [win_left;win_right];

    % Create the sin windows
    case "SIN"
        win_short = sin(pi/N_short * ((0:N_short-1) + 1/2))';
        win_long = sin(pi/N_long * ((0:N_long-1) + 1/2))';
end

% Create the window according to the frame type
switch frameType
    case "OLS"

        % For OLS frame types the window to be applied is a long window
        window = win_long;
    case "ESH"

        % For ESH frame types the short window is repeated 8 times so that
        % it can be applied at the same time in all the subframes
        window = win_short;
        window = repmat(window,8,1);
    case "LSS"

        % If the frame type is LSS then a custom window is created that has
        % 1024 samples from the long window(the first 1024) 448 ones 128
        % (the last 128) from the short window and 448 zeros at the end
        window = [win_long(1:1024);ones(448,1);win_short(129:end);zeros(448,1)];
    case "LPS"

        % If the frame type is LPS then a custom window is created with 448
        % zeros at the end 128 samples from the short window(the first 128)
        % 448 ones and 1024 samples from the long window(the last 1024)
        window = [zeros(448,1);win_short(1:128);ones(448,1);win_long(1025:end)];
end



if (frameType == "ESH")

    % If the frame is eight short then calculate the inverse mdct for every
    % sub frame add the overlapping parts of the frames and then zero pad
    % at the beginning and at the end
    frameF = [reshape(frameF(:,1:2:15),[],1) reshape(frameF(:,2:2:16),[],1)];
    counter_bottom = 1;
    counter_top = counter_bottom + 127;

    frameT = [];
    for i = 1:8
        frameT = [frameT; W_short * frameF(counter_bottom:counter_top,:)];
        counter_bottom = counter_bottom + 128;
        counter_top = counter_top + 128;
    end

    % Apply the window to all the subframes at once
    frameT = frameT .* window;

    channel1 = reshape(frameT(:,1),[256,8]);
    channel2 = reshape(frameT(:,2),[256,8]);
    total1 = zeros(2048,1);
    total2 = zeros(2048,1);
    counter_bottom = 449 + 128;
    total1(449:449+255) = channel1(:,1);
    total2(449:449+255) = channel2(:,1);

    % Add the overlapping parts
    for i = 1:7
        total1(counter_bottom:counter_bottom+127) = total1(counter_bottom:counter_bottom+127) + channel1(1:128,i+1);
        total2(counter_bottom:counter_bottom+127) = total2(counter_bottom:counter_bottom+127) + channel2(1:128,i+1);
        counter_bottom = counter_bottom + 128;
        total1(counter_bottom:counter_bottom+127) = channel1(129:end,i+1);
        total2(counter_bottom:counter_bottom+127) = channel2(129:end,i+1);
    end
    frameT = [total1 total2];
else
    
    frameT = NaN(N_long,2);
    
    % For every other frame just calculate the inverse mdct and multiply
    % the result with the corresponding window
    
    frameT = W_long * frameF;
    frameT = frameT .* window;
end

end