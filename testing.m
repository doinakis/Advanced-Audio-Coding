clc 
clear

[y,Fs] = audioread('./files/LicorDeCalandraca.wav');
frame1 = NaN(2048,ceil(size(y,2)/2048));
frame2 = NaN(2048,ceil(size(y,2)/2048));
counter_top = 2048;
counter_bottom = 1;
for i = 1:2*ceil(size(y,1)/2048)+1
    
    if i == 1
        frame1(:,i) = [zeros(1024,1);y(1:1024,1)];
        frame2(:,i) = [zeros(1024,1);y(1:1024,2)];
        continue;
    end
    frame1(:,i) = y(counter_bottom:counter_top,1);
    frame2(:,i) = y(counter_bottom:counter_top,2);
    counter_top = counter_top + 1024;
    counter_bottom = counter_bottom + 1024;
    if counter_top > size(y,1)
        frame1(:,i) = [y(counter_bottom:end,1);zeros(counter_top - size(y,1),1)];
        frame2(:,i) = [y(counter_bottom:end,2);zeros(counter_top - size(y,1),1)];
        break;
    end
end


b = [0.7548 -0.7548];
a = [1 -0.5095];
frameT = [frame1(:,1) frame2(:,1)];
filtered = filter(b,a,frameT,[],1);

s_values = NaN(8,2);
counter_bottom = 577;
counter_top = 685;
for i = 1:8

    s(i,1) = sum(filtered(counter_bottom:counter_top,1).^2);
    s(i,2) = sum(filtered(counter_bottom:counter_top,2).^2);
    counter_bottom = counter_bottom + 128;
    counter_top = counter_top + 128;
end

attack_values = s./((1/8)* sum(s,1));

counter_bottom = 449;
counter_top = 704;
for i = 1:8
   subframe1(:,i) = frameT(counter_bottom:counter_top,1);
   subframe2(:,i) = frameT(counter_bottom:counter_top,2);
   counter_bottom = counter_bottom + 128;
   counter_top = counter_top + 128;
end
