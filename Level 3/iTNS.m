%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doinakis Michail
% doinakis@ece.auth.gr
% 9292
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that calculates one channel's frame(frameFout given the MDCT
% coeffs after the TNS the frameType and the TNS coefficients.
% Where:
% frameFin: The MDCT coeffs after TNS is applied to them
% frameType: The type of the frame (OLS,LPS,LSS,ESH)
% TNScoeffs: The TNS coefficients
%%
function frameFout = iTNS(frameFin,frameType,TNScoeffs)

if frameType == "ESH"
    % If the frame is ESH then the output frame is 128-by-8
    frameFout = NaN(128,8);
    for i = 1:8
        % Apply the inverse filter
        frameFout(:,i) = filter(1,[1;-TNScoeffs(:,i)],frameFin(:,i),[],1);
    end
else
    % Apply the inverse filter
    frameFout = filter(1,[1;-TNScoeffs],frameFin,[],1);
end
end

