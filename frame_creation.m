
[y,Fs] = audioread('./files/LicorDeCalandraca.wav','double');
frame_counter = 1;
frame1 = NaN(2048,ceil(size(y,1)/1024));
frame2 = NaN(2048,ceil(size(y,1)/1024));
y = [zeros(1024,2);y];

for i = 1:1024:size(y,1)
    if i+2047 > size(y,1)
        frame1(1:2048,frame_counter) = [y(i:end,1);zeros(i+2047-size(y,1),1)];
        frame2(1:2048,frame_counter) = [y(i:end,2);zeros(i+2047-size(y,1),1)];
        break;
    end
    frame1(1:2048,frame_counter) = y(i:i+2047,1);
    frame2(1:2048,frame_counter) = y(i:i+2047,2);
    frame_counter = frame_counter + 1;
end

frameTypes = strings(frame_counter,1);
for i = 1:frame_counter
    frameT = [frame1(:,i) frame2(:,i)];
    
    if i == frame_counter 
        nextframeT = zeros(2048,2);
    else
        nextframeT = [frame1(:,i+1) frame2(:,i+1)];
    end
    if i == 1
        prevframeType = "OLS";
        frameTypes(i) = SSC(frameT,nextframeT,prevframeType);
        continue;
    end
    prevframeType = frameTypes(i-1);
    frameTypes(i) = SSC(frameT,nextframeT,prevframeType);
end

clearvars  frame_counter i prevframeType frameT frame_counter y;