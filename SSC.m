%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OLS = ONLY_LONG_SEQUENCE
%   LSS = LONG_START_SEQUENCE
%   ESH = EIGHT_SHORT_SEQUENCE
%   LPS = LONG_STOP_SEQUENCE
%%


function frameType = SSC(frameT,nextframeT,prevframeType)

H = @(z) (0.7548 - 0.7548 .* z^(-1))/(1 - 0.5095 .* z^(-1));

switch prevframeType
    case "LSS"
        frameType = "ESH";
        return;
    case "LPS"
        frameType = "OLS";
        return;
    case "OLS"
        
    case "ESH"
end


end

