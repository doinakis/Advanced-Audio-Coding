clc 
clear

[y,Fs] = audioread('./files/LicorDeCalandraca.wav','double');
frame_counter = 1;
for i = 1:1024:size(y,1)
%     if i == 1 
%         frame1(:,frame_counter)=[zeros(1024,1);y(1:1024,1)];
%         frame2(:,frame_counter)=[zeros(1024,1);y(1:1024,2)];
%         frame_counter = frame_counter + 1;
%         continue;
%     end
    if i+2047 > size(y,1)
        frame1(1:2048,frame_counter) = [y(i:end,1);zeros(i+2047-size(y,1),1)];
        frame2(1:2048,frame_counter) = [y(i:end,2);zeros(i+2047-size(y,1),1)];
        break;
    end
    frame1(1:2048,frame_counter) = y(i:i+2047,1);
    frame2(1:2048,frame_counter) = y(i:i+2047,2);
    frame_counter = frame_counter + 1;
end

frameTypes = strings(276,1);
for i = 1:276
    frameT = [frame1(:,i) frame2(:,i)];
    
    if i == 276 
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
signal = [];
for i = 1:276
    signal =[signal;frame1(:,i) frame2(:,i)];
end