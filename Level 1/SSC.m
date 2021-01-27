%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that calculates the i-th frames type based on the previous frame
% and the next frame. It returns the frame type as follows
% "OLS" = ONLY_LONG_SEQUENCE
% "LSS" = LONG_START_SEQUENCE
% "ESH" = EIGHT_SHORT_SEQUENCE
% "LPS" = LONG_STOP_SEQUENCE
% Where:
% frameT: 2048-by-2 array with the values of the current frame
% nextframeT: 2048-by-2 array with the values of the next frame
% prevframeType: The type of the previous frame
%%
function frameType = SSC(frameT,nextframeT,prevframeType)

% Initialize the filter parameters
a = [1 -0.5095];
b = [0.7548 -0.7548];

% Find the frame Type acording to the algorithm specified
switch prevframeType

    % If the previous frame is LSS the next frame is ESH
    case "LSS"
        frameType = "ESH";
        return;

    % If the previous frame is LPS the next frame is OLS
    case "LPS"
        frameType = "OLS";
        return;
    % If the previous frame is OLS the next frame should me studied first
    case "OLS"

        % Apply the filter to the next frame to check if it is ESH
        % The filter applied: H(z) = (0.7548
        % -7548*z^(-1))/(1-0.5095*z^(-1))
        filtered_samples = filter(b,a,nextframeT,[],1);

        % Preallocate space for the s and attack values
        attack_values = NaN(7,2);
        s_values = NaN(8,2);
        counter = 1;

        % The bottom counter starts from the 576-th sample of the nextframe
        counter_bottom = 449 + 128;
        counter_top = counter_bottom + 127;

        % Initial values for the left and right channel
        frameType1 = "OLS";
        frameType2 = "OLS";

        % Find the frame type of each channel by calculating the s and
        % attack values
        for i = 1:8
            s_values(i,:) = sum(filtered_samples(counter_bottom:counter_top,:).^2);
            if i > 1
                attack_values(i-1,:) = (counter) * s_values(i,:)./sum(s_values(1:counter,:),1);
                counter = counter + 1;

                % Check for both channelsthe attack and s values to see if
                % the next frame is ESH or not if it is then the current
                % channels frame is LSS else the initial values remains
                if (s_values(i,1) > 1e-3 && attack_values(i-1,1) > 10)
                    frameType1 = "LSS";
                end
                if (s_values(i,2) > 1e-3 && attack_values(i-1,2) > 10)
                    frameType2 = "LSS";
                end
            end

            % Increase the counters to check the next subframe
            counter_bottom = counter_bottom + 128;
            counter_top = counter_top + 128;
        end
    case "ESH"

        % If the previous frame is ESH the next frame should me studied first
        % The filter applied: H(z) = (0.7548
        % -7548*z^(-1))/(1-0.5095*z^(-1))
        filtered_samples = filter(b,a,nextframeT,[],1);

        % Preallocate space for the s and attack values
        attack_values = NaN(7,2);
        s_values = NaN(8,2);
        counter = 1;
        counter_bottom = 449 + 128;
        counter_top = counter_bottom + 127;

        % Initial values for the left and right channel
        frameType1 = "LPS";
        frameType2 = "LPS";

        % Calculations for every subframe
        for i = 1:8
            s_values(i,:) = sum(filtered_samples(counter_bottom:counter_top,:).^2);
            if i > 1
                attack_values(i-1,:) = (counter) * s_values(i,:)./sum(s_values(1:counter,:),1);
                counter = counter + 1;

                % Check for both channelsthe attack and s values to see if
                % the next frame is ESH or not if it is then the current
                % channels frame is ESH else the initial values remains
                if (s_values(i,1) > 1e-3 && attack_values(i-1,1) > 10)
                    frameType1 = "ESH";
                end
                if (s_values(i,2) > 1e-3 && attack_values(i-1,2) > 10)
                    frameType2 = "ESH";
                end
            end

            % Increase the counters to check the next subframe
            counter_bottom = counter_bottom + 128;
            counter_top = counter_top + 128;
        end
end

% Specify the final frameType using the frame types of each channel
% calculated above
% If the left and right channel types match each other the the frame is the
% same as those
if frameType1 == frameType2
    frameType = frameType1;
    return;
end

% If one of them is OLS and the other LSS then the frame type is LSS
if (frameType1 == "OLS" && frameType2 == "LSS") || (frameType1 == "LSS" && frameType2 == "OLS")
    frameType = "LSS";
    return;

% If one of them is OLS and the other LPS then the frame type is LPS
elseif (frameType1 == "OLS" && frameType2 == "LPS") || (frameType1 == "LPS" && frameType2 == "OLS")
    frameType = "LPS";
    return;

% In every other case the frame type is ESH
else
    frameType = "ESH";
    return;
end
end

