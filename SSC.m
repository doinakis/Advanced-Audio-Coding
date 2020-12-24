%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OLS = ONLY_LONG_SEQUENCE
%   LSS = LONG_START_SEQUENCE
%   ESH = EIGHT_SHORT_SEQUENCE
%   LPS = LONG_STOP_SEQUENCE
%%
function frameType = SSC(frameT,nextframeT,prevframeType)


a = [1 -0.5095];
b = [0.7548 -0.7548];

switch prevframeType
    case "LSS"
        frameType = "ESH";
        return;
    case "LPS"
        frameType = "OLS";
        return;
    case "OLS"
        filtered = filter(b,a,nextframeT,[],1);
        s_values = NaN(8,2);
        counter_bottom = 577;
        counter_top = 685;
        for i = 1:8
            s_values(i,1) = sum(filtered(counter_bottom:counter_top,1).^2);
            s_values(i,2) = sum(filtered(counter_bottom:counter_top,2).^2);
            counter_bottom = counter_bottom + 128;
            counter_top = counter_top + 128;
        end
        attack_values = s_values./((1/8)* sum(s_values,1));
        
        if any(s_values(:,1) > 1e-3) && any(attack_values(:,1) > 10)
            frame_channel_1 = "LSS";
        else
            frame_channel_1 = "OLS";
        end
        if any(s_values(:,2) > 1e-3) && any(attack_values(:,2) > 10)
            frame_channel_2 = "LSS";
        else
            frame_channel_2 = "OLS";
        end
        if frame_channel_1 == frame_channel_2 
            frameType = frame_channel_1;
        elseif (frame_channel_1 == "OLS" && frame_channel_2 == "LSS") || (frame_channel_2 == "OLS" && frame_channel_1 == "LSS")
            frameType = "LSS";
        elseif (frame_channel_1 == "OLS" && frame_channel_2 == "LPS") || (frame_channel_2 == "OLS" && frame_channel_1 == "LPS")
            frameType = "LPS";
        else
            frameType = "ESH";
        end
    case "ESH"
        filtered = filter(b,a,nextframeT,[],1);
        s_values = NaN(8,2);
        counter_bottom = 577;
        counter_top = 685;
        for i = 1:8
            s_values(i,1) = sum(filtered(counter_bottom:counter_top,1).^2);
            s_values(i,2) = sum(filtered(counter_bottom:counter_top,2).^2);
            counter_bottom = counter_bottom + 128;
            counter_top = counter_top + 128;
        end
        attack_values = s_values./((1/8)* sum(s_values,1));
        
        if any(s_values(:,1) > 1e-3) && any(attack_values(:,1) > 10)
            frame_channel_1 = "ESH";
        else
            frame_channel_1 = "LPS";
        end
        if any(s_values(:,2) > 1e-3) && any(attack_values(:,2) > 10)
            frame_channel_2 = "ESH";
        else
            frame_channel_2 = "LPS";
        end
        if frame_channel_1 == frame_channel_2 
            frameType = frame_channel_1;
        elseif (frame_channel_1 == "OLS" && frame_channel_2 == "LSS") || (frame_channel_2 == "OLS" && frame_channel_1 == "LSS")
            frameType = "LSS";
        elseif (frame_channel_1 == "OLS" && frame_channel_2 == "LPS") || (frame_channel_2 == "OLS" && frame_channel_1 == "LPS")
            frameType = "LPS";
        else
            frameType = "ESH";
        end
end
end

