function frameFout = iTNS(frameFin,frameType,TNScoeffs)

if frameType == "ESH"
    frameFout = NaN(128,8);
    for i = 1:8
        frameFout(:,i) = filter(1,[1;TNScoeffs(:,i)],frameFin(:,i),[],1);
    end
else
    frameFout = filter(1,[1;TNScoeffs],frameFin,[],1);
end
end

