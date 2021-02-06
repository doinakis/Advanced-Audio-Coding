%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that implements the filterbank. It returns the mdct coefficients
% of frameT. If the frame is ESH then the frameT contains 8 128-by-2
% matrices (added in columns), otherwise is contains a 1024-by-2.
% Where:
% frameT: The current frame to calculate its coefficients
% frameType: The type of the frame (OLS,LPS,LSS,ESH)
% winType: The type of the window to be applied (SIN,KBD)
%%
function frameF = filterbank(frameT,frameType,winType)

% MDCT matrices coefficients
persistent W_long W_short

if isempty(W_long)
    % Calculate the MDCT coefficients matrices for short and long windows
    N_long = 2048;
    n_zero = (N_long/2 + 1)/2;
    n = 0:N_long-1;
    k = (0:N_long/2-1)';
    W_long = 2*cos(2*pi/N_long .* (n + n_zero) .* (k + 1/2));
    
    N_short = 256;
    n_zero = (N_short/2 + 1)/2;
    n = 0:N_short-1;
    k = (0:N_short/2-1)';
    W_short = 2*cos(2*pi/N_short .* (n + n_zero) .* (k + 1/2));
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
    
    % If the frame is ESH then calculate the 8 areas
    y = NaN(2048,2);
    
    % Initialize an array for the mdct coefficients
    frameF = [];
    % The y variable holds 2048-by-2 samples. It holds the 8*256 subframes
    counter = 1;
    for i = 449:128:1345
        y(counter:counter+255,:) = frameT(i:i+255,:);
        counter = counter + 256;
    end

    % Apply the window to every subframe at once
    z = y .* window;

    % Apply the mdct
    counter_bottom = 1;
    counter_top = counter_bottom + 255;
    
    for i = 1:8
        frameF = [frameF W_short * z(counter_bottom:counter_top,:)];
        counter_bottom = counter_bottom + 256;
        counter_top = counter_top + 256;
    end

else
    
    % Initialize an array for the mdct coefficients
    frameF = NaN(N_long/2,2);
    
    % Apply the window to the current frame
    z = frameT .* window;

    % Calculate the mdct
    frameF = W_long * z;
end

end

