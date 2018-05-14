function Project_2()
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
Nlp = zeros(128, 2);
Dlp = zeros(128, 2);
Np = zeros(128, 3);
Dp = zeros(128, 3);
Nz = zeros(128, 3);
Dz = zeros(128, 3);

for i=1:128
    for j=1:3
        if j == 1
            Nlp(i,j) = 1-a0(i);
            Dlp(i,j) = 1;
            Np(i,j) = 1-b1(i)+b2(i) ;       
            Dp(i,j) = 1;
            Nz(i,j) = 1;        
            Dz(i,j) = 1-a1(i)+a2(i);  
        else if j == 2
            Dlp(i,j) = -a0(i);
            Dp(i,j) = -b1(i);
            Nz(i,j) = -a1(i);        
            else
            Dp(i,j) = b2(i);
            Nz(i,j) = a2(i) ;       
            end
        end
    end
end

% % low pass filter frequency response
% figure;
% for k=1:16
%     filter_no = k*8;    
%     freqz([Nlp(filter_no,1) Nlp(filter_no,2)],[Dlp(filter_no,1) Dlp(filter_no,2)], F, fs);
%     ax = findall(gcf, 'Type', 'axes');
%     set(ax, 'XScale', 'log');
%     hold on
% end
% 
% resonant pole frequency response
% figure;
% for k=1:16
%     filter_no = k*8;
%     freqz([Np(filter_no,1) Np(filter_no,2) Np(filter_no,3)],[Dp(filter_no,1) Dp(filter_no,2) Dp(filter_no,3)], F, fs);
%     ax = findall(gcf, 'Type', 'axes');
%     set(ax, 'XScale', 'log');
%     hold on
% end
% 
% % resonant zero frequency response
% figure;
% for k=1:16
%     filter_no = k*8;    
%     freqz([Nz(filter_no,1) Nz(filter_no,2) Nz(filter_no,3)],[Dz(filter_no,1) Dz(filter_no,2) Dz(filter_no,3)], F, fs);
%     ax = findall(gcf, 'Type', 'axes');
%     set(ax, 'XScale', 'log');
%     hold on
% end

%% displacement transfer function
RLP_num = zeros(128, 4);
RLP_den = zeros(128, 4);

% cascade of low pass and resonant pole transfer functions
for i=1:128
    RLP_num(i,:) = conv([Nlp(i,1) Nlp(i,2)],[Np(i,1) Np(i,2) Np(i,3)]);
    RLP_den(i,:) = conv([Dlp(i,1) Dlp(i,2)],[Dp(i,1) Dp(i,2) Dp(i,3)]);
end

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


%% Filters
figure;
samples = 1:20000;
x1 = [1, zeros(1, length(samples) - 1)];
x2 = sin(2*pi*1000*samples/fs);
x3 = sin(2*pi*16000*samples/fs) + sin(2*pi*4000*samples/fs) + sin(2*pi*200*samples/fs) ;
t = linspace(0, 1000*length(samples)/fs, length(samples));

for filter_no = 1:128;
    if filter_no == 1
        d1(filter_no,:) = filter([RLP_num(filter_no,:)],[RLP_den(filter_no,:)], x1);
        d2(filter_no,:) = filter([RLP_num(filter_no,:)],[RLP_den(filter_no,:)], x2);
        d3(filter_no,:) = filter([RLP_num(filter_no,:)],[RLP_den(filter_no,:)], x3);
        p1(filter_no,:) = filter([RNLP_num(filter_no,:)],[RNLP_den(filter_no,:)], x1);
        p2(filter_no,:) = filter([RNLP_num(filter_no,:)],[RNLP_den(filter_no,:)], x2);
        p3(filter_no,:) = filter([RNLP_num(filter_no,:)],[RNLP_den(filter_no,:)], x3);
    else
        d1(filter_no,:) = filter([RLP_num(filter_no,:)],[RLP_den(filter_no,:)], d1(filter_no - 1,:));
        d2(filter_no,:) = filter([RLP_num(filter_no,:)],[RLP_den(filter_no,:)], d2(filter_no - 1,:));
        d3(filter_no,:) = filter([RLP_num(filter_no,:)],[RLP_den(filter_no,:)], d3(filter_no - 1,:));
        p1(filter_no,:) = filter([RNLP_num(filter_no,:)],[RNLP_den(filter_no,:)], p1(filter_no - 1,:));
        p2(filter_no,:) = filter([RNLP_num(filter_no,:)],[RNLP_den(filter_no,:)], p2(filter_no - 1,:));
        p3(filter_no,:) = filter([RNLP_num(filter_no,:)],[RNLP_den(filter_no,:)], p3(filter_no - 1,:));

    end
end


%% Spatial differentiation
delta_x = 3.5/128;

for i = 1:128
    if i < 128
    s1(i,:) = (d1(i+1,:) - d1(i,:)) / delta_x;
    s2(i,:) = (d2(i+1,:) - d2(i,:)) / delta_x;
    s3(i,:) = (d3(i+1,:) - d3(i,:)) / delta_x;
    else
        s1(i,:) = d1(i,:);
        s2(i,:) = d2(i,:);
        s3(i,:) = d3(i,:);
    end
end
s(1,:) = s2(1,:);
for i=2:128    
    s(i,:) = s(i-1,:) +s2(i,:);
    % title('Impulse response of filter');
    % xlabel('Time (ms)');
    %ylabel('Displacement');
end
plot(s(81,:))
end