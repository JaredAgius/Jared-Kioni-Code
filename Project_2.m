function Project_2(Xin)
format compact
close all
clc

fs = 48000;
x = linspace(0, 3.5, 128);  %distance
z = 1.0432; 
F = 10:20000;

qp = linspace(10, 5.5, 128); %pole quality factor
qz = linspace(22, 12, 128); %zeros quality factor

%% resonant freq, notch filter freq, bandwidth (z & p)
for i = 1:128
    fp(i) = 20000 * 10 .^(-0.667*x(i));
    fz(i) = fp(i) * z;
    bwp(i) = fp(i)./qp(i);
    bwz(i) = fz(i)./qz(i);
    tp(i) = 2 * pi * fp(i) / fs;
    tz(i) = 2 * pi * fz(i) / fs;
    rp(i) = 1 - (bwp(i)/fs) * pi;
    rz(i) = 1 - (bwz(i)/fs) * pi;
    b1(i) = 2 * rp(i) * cos(tp(i));
    b2(i) = rp(i).^2;
    a1(i) = 2 * rz(i) * cos(tz(i));
    a2(i) = rz(i).^2;
    flp(i) = fz(i) * 1.4;
    tlp(i) = 2 * pi * flp(i) / fs;
    a0(i) = (2 - cos(tlp(i))) - sqrt((2-cos(tlp(i))).^2 - 1);
end

%% Low pass filter, resonant pole, resonant zero transfer functions
Nlp = zeros(128, 2);    % Low Pass Numerator
Dlp = zeros(128, 2);    % Low Pass Denominator
Np = zeros(128, 3);     % Resonant Pole Numerator
Dp = zeros(128, 3);     % Resonant Pole Denominator
Nz = zeros(128, 3);     % Resonant Zero Numerator
Dz = zeros(128, 3);     % Resonant Zero Denominator
Npd = zeros(128, 3);    % Resonant Zero Numerator for Displacement output

for i=1:128
    for j=1:3
        if j == 1
            Nlp(i,j) = 1-a0(i);
            Dlp(i,j) = 1;
            Np(i,j) = 1-b1(i)+b2(i);    
            Dp(i,j) = 1;
            Nz(i,j) = 1;        
            Dz(i,j) = 1-a1(i)+a2(i);  
        else if j == 2
            Dlp(i,j) = -a0(i);
            Dp(i,j) = -b1(i);
            Nz(i,j) = -a1(i);
            Npd(i,j) = 1-b1(i)+b2(i);
        else
            Dp(i,j) = b2(i);
            Nz(i,j) = a2(i) ;       
            end
        end
    end
end

%% Low pass filter, resonant pole, resonant zero frequency response
%{
% Low pass filter frequency response
 figure;
 for k=1:16
     filter_no = k*8;    
     freqz([Nlp(filter_no,1) Nlp(filter_no,2)],[Dlp(filter_no,1) Dlp(filter_no,2)], F, fs);
     ax = findall(gcf, 'Type', 'axes');
     set(ax, 'XScale', 'log');
     hold on
 end
 
 % Resonant pole frequency response
 figure;
 for k=1:16
     filter_no = k*8;
     freqz([Np(filter_no,1) Np(filter_no,2) Np(filter_no,3)],[Dp(filter_no,1) Dp(filter_no,2) Dp(filter_no,3)], F, fs);
     ax = findall(gcf, 'Type', 'axes');
     set(ax, 'XScale', 'log');
     hold on
 end

% Resonant zero frequency response
 figure;
 for k=1:16
     filter_no = k*8;    
     freqz([Nz(filter_no,1) Nz(filter_no,2) Nz(filter_no,3)],[Dz(filter_no,1) Dz(filter_no,2) Dz(filter_no,3)], F, fs);
     ax = findall(gcf, 'Type', 'axes');
     set(ax, 'XScale', 'log');
     hold on
 end
%}

%% displacement transfer function


RLP_num = zeros(128, 4);
RLP_den = zeros(128, 4);

% cascade of low pass and resonant pole transfer functions
for i=1:128
    RLP_num(i,:) = conv([Nlp(i,1) Nlp(i,2)],[Np(i,2) Np(i,1) Np(i,3)]);     % SWAP 1 and 2
    RLP_den(i,:) = conv([Dlp(i,1) Dlp(i,2)],[Dp(i,1) Dp(i,2) Dp(i,3)]);
end

figure
freqz(RLP_num(81,:),RLP_den(81,:))
title('Displacement TF')
ax = findall(gcf, 'Type', 'axes');
set(ax, 'XScale', 'log');

%{
% figure;
% %displacement transfer function frequency resposne
% for k=1:16
%     filter_no = k*8;
%     freqz([RLP_num(filter_no,1) RLP_num(filter_no,2) RLP_num(filter_no,3) RLP_num(filter_no,4)],[RLP_den(filter_no,1) RLP_den(filter_no,2) RLP_den(filter_no,3) RLP_den(filter_no,4)], F, fs);
%     ax = findall(gcf, 'Type', 'axes');
%     set(ax, 'XScale', 'log');
%     hold on
% end

%% pressure transfer function
RNLP_num = zeros(128, 6);
RNLP_den = zeros(128, 6);

%cascade of low pass filter, resonant pole filter and resonant zero filter
for i=1:128
    RNLP_num(i,:) = conv([RLP_num(i,1) RLP_num(i,2) RLP_num(i,3) RLP_num(i,4)],[Nz(i,1) Nz(i,2) Nz(i,3)]);
    RNLP_den(i,:) = conv([RLP_den(i,1) RLP_den(i,2) RLP_den(i,3) RLP_den(i,4)],[Dz(i,1) Dz(i,2) Dz(i,3)]);
end
 
% % pressure transfer function frequency response
% figure;
% for k=1:16
%     filter_no = k*8;
%     freqz([RNLP_num(filter_no,1) RNLP_num(filter_no,2) RNLP_num(filter_no,3) RNLP_num(filter_no,4) RNLP_num(filter_no,5) RNLP_num(filter_no,6)],[RNLP_den(filter_no,1) RNLP_den(filter_no,2) RNLP_den(filter_no,3) RNLP_den(filter_no,4) RNLP_den(filter_no,5) RNLP_den(filter_no,6)], F, fs);
%     ax = findall(gcf, 'Type', 'axes');
%     set(ax, 'XScale', 'log');
%     hold on
% end

%}
%% Filters
samples = 1:8000;
t = samples/fs;

% Inputs (Dependant on Xin)
x1 = zeros(3, length(samples));
x1(1,:) = [1, zeros(1, length(samples) - 1)];
x1(2,:) = sin(2*pi*1000*samples/fs);
x1(3,:) = sin(2*pi*16000*samples/fs) + sin(2*pi*4000*samples/fs) + sin(2*pi*200*samples/fs) ;

for filter_no = 1:128;
   if filter_no == 1
        ylp(filter_no,:) = filter([Nlp(filter_no,:)],[Dlp(filter_no,:)],x1(Xin,:));
        yrp(filter_no,:) = filter([Np(filter_no,:)],[Dp(filter_no,:)],ylp(filter_no,:));
        yrz(filter_no,:) = filter([Nz(filter_no,:)],[Dz(filter_no,:)],yrp(filter_no,:));
        ydisp(filter_no,:) = filter([Npd(filter_no,:)],[Dp(filter_no,:)],ylp(filter_no,:));
   else
        ylp(filter_no,:) = filter([Nlp(filter_no,:)],[Dlp(filter_no,:)],yrz(filter_no-1,:));
        yrp(filter_no,:) = filter([Np(filter_no,:)],[Dp(filter_no,:)],ylp(filter_no,:));
        yrz(filter_no,:) = filter([Nz(filter_no,:)],[Dz(filter_no,:)],yrp(filter_no,:));
        ydisp(filter_no,:) = filter([Npd(filter_no,:)],[Dp(filter_no,:)],ylp(filter_no,:));
   end    
end

% Pressure Envelope
plot(yrz);xlabel('Samples');ylabel('Disp');

% Impulse Response 
figure
f(1) = subplot(3,1,1);
plot(t*1000,ydisp(30,:))
f(2) = subplot(3,1,2);
plot(t*1000,ydisp(60,:))
f(3) = subplot(3,1,3);
plot(t*1000,ydisp(90,:))

title(f(1),'Impulse Response Filter 30');ylabel(f(1),'Disp.');xlabel(f(1),'Time (ms)')
title(f(2),'Impulse Response Filter 60');ylabel(f(2),'Disp.');xlabel(f(2),'Time (ms)')
title(f(3),'Impulse Response Filter 90');ylabel(f(3),'Disp.');xlabel(f(3),'Time (ms)')

% Frequency Response
figure
nfft = 8000;
YDISP30 = dtft(ydisp(30,:),nfft);
YDISP60 = dtft(ydisp(60,:),nfft);
YDISP90 = dtft(ydisp(90,:),nfft);
f(1) = subplot(3,1,1);
plot(abs(YDISP30))
f(2) = subplot(3,1,2);
plot(abs(YDISP60))
f(3) = subplot(3,1,3);
plot(abs(YDISP90))

title(f(1),'Frequency Response Filter 30');ylabel(f(1),'Disp. Ratio');xlabel(f(1),'Samples')
title(f(2),'Frequency Response Filter 60');ylabel(f(2),'Disp. Ratio');xlabel(f(2),'Samples')
title(f(3),'Frequency Response Filter 90');ylabel(f(3),'Disp. Ratio');xlabel(f(3),'Samples')

%% Spatial differentiation
delta_x = 3.5/128;

for i = 1:128
    if i < 128
        s1(i,:) = (ydisp(i+1,:) - ydisp(i,:)) / delta_x;
    else
        s1(i,:) = ydisp(i,:);
    end
end

for i = 1:128
    if i < 128
        s2(i,:) = (s1(i+1,:) - s1(i,:)) / delta_x;
    else
        s2(i,:) = s1(i,:);
    end
end

% Effect of Spatial Differentiation
 DispIndex = 384;
 figure;
 f(1) = subplot(3,1,1);
 plot(ydisp(:,DispIndex))
 f(2) = subplot(3,1,2);
 plot(s1(:,DispIndex))
 f(3) = subplot(3,1,3);
 plot(s2(:,DispIndex))
 
 title(f(1),'Displacement Output: 0 Spacial Differnetiation');ylabel(f(1),'Disp.');xlabel(f(1),'Filter no.')
 title(f(2),'Displacement Output: 1 Spacial Differnetiation');ylabel(f(1),'Disp.');xlabel(f(1),'Filter no.')
 title(f(3),'Displacement Output: 2 Spacial Differnetiations');ylabel(f(1),'Disp.');xlabel(f(1),'Filter no.')
 
 %% Inner Hair Cell Model
 Co = exp(-30*2*pi/48000);
 
 % Half wave rectification 
 for i = 1:128
    if s2(i,DispIndex) < 0
       s2(i,DispIndex) = 0; 
    end
 end

 % Output Electrical Signal
 for i = 1:128
    if i == 1
        v(i) = (1 - Co)*s2(i,DispIndex);
    else
        v(i) = (1 - Co)*s2(i,DispIndex) + Co*v(i-1);
    end
 end

 figure
 plot(v);ylabel('Disp.');xlabel('Filter no.')
 
 
end
