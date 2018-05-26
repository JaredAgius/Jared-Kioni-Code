% ELEC3104: Project - Cochlear Signal Processing
% Authors:  Jared Agius     z5113699
%           Kioni Ndirangu  z5111826
% Date:     27/05/18

% This code displays a MATLAB implementation of the transmission line cochlea model as digital filters. 
% The response of the basillar membrane to differnt inputs will be displayed and the implementation will be validated in terms of the 
% systems impulse and magnitude response

function Project_2(Input,OuterMiddle)
format compact
close all
clc

fs = 48000;                                                         % Sampling frequency 
x = linspace(0, 3.5, 128);                                          % Distance of filters
z = 1.0432;                                                         % Pole - Zero frequency relation           
qp = linspace(10, 5.5, 128);                                        % Poles quality factor
qz = linspace(22, 12, 128);                                         % Zeros quality factor

%% Resonant, Notch and Low Pass Filter Frequncies and Coefficients. Bandwidth of Poles and Zeros

for i = 1:128
    fp(i) = 20000 * 10 .^(-0.667*x(i));                             % Pole Frequencis
    fz(i) = fp(i) * z;                                              % Zero Frequencies
    bwp(i) = fp(i)./qp(i);                                          % Pole Bandwiths
    bwz(i) = fz(i)./qz(i);                                          % Zero Bandwiths
    tp(i) = 2 * pi * fp(i) / fs;                                    % Pole Digital Frequncies
    tz(i) = 2 * pi * fz(i) / fs;                                    % Zero Digital Frequencies
    rp(i) = 1 - (bwp(i)/fs) * pi;                                   % Pole Radius
    rz(i) = 1 - (bwz(i)/fs) * pi;                                   % Zero Radius
    b1(i) = 2 * rp(i) * cos(tp(i));                                 % B1 Coeffiecint
    b2(i) = rp(i).^2;                                               % B2 Coeffiecint
    a1(i) = 2 * rz(i) * cos(tz(i));                                 % A1 Coeffiecint
    a2(i) = rz(i).^2;                                               % A2 Coeffiecint
    flp(i) = fz(i) * 1.4;                                           % Low Pass Cutoff Frequency 
    tlp(i) = 2 * pi * flp(i) / fs;                                  % Low Pass Digital Frequencies
    a0(i) = (2 - cos(tlp(i))) - sqrt((2-cos(tlp(i))).^2 - 1);       % A0 Coeffiecint
end

%% Low pass filter, resonant pole, resonant zero transfer functions

Nlp = zeros(128, 2);                                                % Low Pass Numerator
Dlp = zeros(128, 2);                                                % Low Pass Denominator
Np = zeros(128, 3);                                                 % Resonant Pole Numerator
Dp = zeros(128, 3);                                                 % Resonant Pole Denominator
Nz = zeros(128, 3);                                                 % Resonant Zero Numerator
Dz = zeros(128, 3);                                                 % Resonant Zero Denominator
Npd = zeros(128, 3);                                                % Resonant Zero Numerator for Displacement output

% Calculate pole, zero and low pass polynomial coeffiecients  
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

%% Cascading Digital Filters
% Inputs to the system pass through a cascade of filters

samples = 1:8000;                                                   % 8000 samples
t = samples/fs;                                                     % Time

% Inputs
x1 = zeros(3, length(samples));
x1(1,:) = [1, zeros(1, length(samples) - 1)];                       % Impulse
x1(2,:) = sin(2*pi*800*samples/fs);                                 % Sinusoid
x1(3,:) = sin(2*pi*16000*samples/fs) + sin(2*pi*4000*samples/fs)... % Multiple Frequency Sinusoid
        + sin(2*pi*200*samples/fs) ;
x1(4,:) = transpose(audioread('BeCool.wav'));

if OuterMiddle                                                      % Pass the input through the outer and middle ear
   x1(Input,:) = MiddleEarGain(x1(Input,:));                        % if OuterMiddle == 1
end    

% Cascaded Filters
% ylp - low pass filter output     
% yrp - Resonant pole filter output
% yrz - Resonant zero filter output (Pressure outputs)
% ydisp - Displacement output

K1 = 0.981;                                                         % Constant to give 0db gain at DC
for filter_no = 1:128;
   if filter_no == 1
        ylp(filter_no,:) = filter([K1*Nlp(filter_no,:)],[Dlp(filter_no,:)],x1(Input,:));
        yrp(filter_no,:) = filter([Np(filter_no,:)],[Dp(filter_no,:)],ylp(filter_no,:));
        yrz(filter_no,:) = filter([Nz(filter_no,:)],[Dz(filter_no,:)],yrp(filter_no,:));
        ydisp(filter_no,:) = filter([Npd(filter_no,:)],[Dp(filter_no,:)],ylp(filter_no,:));
   else
        ylp(filter_no,:) = filter([K1*Nlp(filter_no,:)],[Dlp(filter_no,:)],yrz(filter_no-1,:));
        yrp(filter_no,:) = filter([Np(filter_no,:)],[Dp(filter_no,:)],ylp(filter_no,:));
        yrz(filter_no,:) = filter([Nz(filter_no,:)],[Dz(filter_no,:)],yrp(filter_no,:));
        ydisp(filter_no,:) = filter([Npd(filter_no,:)],[Dp(filter_no,:)],ylp(filter_no,:));
   end    
end

% Pressure Envelope
plot(yrz(:,(4000:end)));xlabel('Samples');ylabel('Disp');

% Displacement Output (Impulse Response for impulse input)
figure
f(1) = subplot(3,1,1);
plot(t*1000,ydisp(30,:))
f(2) = subplot(3,1,2);
plot(t*1000,ydisp(60,:))
f(3) = subplot(3,1,3);
plot(t*1000,ydisp(90,:))

title(f(1),'y[n] Filter 30');ylabel(f(1),'Disp.');xlabel(f(1),'Time (ms)');axis(f(1),[0,25,-0.1,0.1]);
title(f(2),'y[n] Filter 60');ylabel(f(2),'Disp.');xlabel(f(2),'Time (ms)');axis(f(2),[0,50,-0.02,0.02]);
title(f(3),'y[n] Filter 90');ylabel(f(3),'Disp.');xlabel(f(3),'Time (ms)');axis(f(3),[0,100,-0.004,0.004]);

% Displacement Output Frequency Response
figure
nfft = 8000;
YDISP30 = dtft(ydisp(30,:),nfft);
YDISP30 = YDISP30(floor(length(YDISP30)/2):end);
semilogx(mag2db(abs(YDISP30)));

title('y({\theta}) Filter 30');ylabel('Disp. Ratio');xlabel('Freq.')

%% Spatial differentiation
% Spacial differentiation of displacement outputs, used to model fluid coupling of the outputs. 
% This has a high pass effect.

delta_x = 3.5/128;                                                  % Distance between filters

% 1st Spacial Differentiation 
for i = 1:128
    if i < 128
        s1(i,:) = (ydisp(i+1,:) - ydisp(i,:)) / delta_x;
    else
        s1(i,:) = ydisp(i,:);
    end
end

% 2nd Spacial Differentiation 
for i = 1:128
    if i < 128
        s2(i,:) = (s1(i+1,:) - s1(i,:)) / delta_x;
    else
        s2(i,:) = s1(i,:);
    end
end

% The Effect of Spatial Differentiation
DispIndex = 4000;                                                  % High sample number: when output has stabilized
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

figure
nfft = 8000;
K2 = 0.7657^30;
YSD30 = dtft(K2*s2(30,:),nfft);
YSD30 = YSD30(floor(length(YSD30)/2):end);
hold on
semilogx(mag2db(abs(YSD30)),'r');
semilogx(mag2db(abs(YDISP30)),'b');
hold off

title('Y({\theta}) Filter 30');ylabel('Disp. Ratio');xlabel('Freq.')
 
 %% Inner Hair Cell Model
 % Differentiated and rectified displacement is translated into electrical energy by the neural transduction mechanism
 
 % Half wave rectification 
 % Setting all negative values to zero
 for i = 1:128
    for j = 1:length(samples)
        if s2(i,j) < 0
           s2(i,j) = 0; 
        end
    end
 end

 % Output Electrical Signal
 % Differentiated displacement is translated into electrical energy by the neural transduction mechanism
 v = zeros(128,length(samples));
 Co = exp(-30*2*pi/48000);                                          % Hair Cell Model Constant
 for i = 1:128
     for j = 1:length(samples)
        if j == 1
            v(i,j) = (1 - Co)*s2(i,j);
        else
            v(i,j) = (1 - Co)*s2(i,j) + Co*v(i,j-1);
        end
    end
 end
 
% Display Inner Hair Cell output
figure
plot(v(:,DispIndex));ylabel('Disp.');xlabel('Filter no.')
 
% End of main program
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELEC3104: Project - Cochlear Signal Processing
% Authors:  Jared Agius     z5113699
%           Kioni Ndirangu  z5111826
% Date:     27/05/18

% This code implements the transfer functions for the outer and middle ear,
% by placing poles and zeros

function [Output] = MiddleEarGain(Input)
%% Middle Ear
%Zeros
r1 = 0.55;
r2 = 0.92;
o1 = 2*pi*20/48;     
o2 = 2*pi*0.1/48; 

B1 = 4*conv([1,-2*r1*cos(o1),r1^2],[1,-2*r2*cos(o2),r2^2]);

%Poles
rp = 0.9;
op = 2*pi*1/48;

A1 = [1 -rp*2*cos(op) rp^2];

%% Outer Ear
%Zeros
ro1 = 0.9;
ro2 = 0.9;
oo1 = 2*pi*0.1/48;           
oo2 = 2*pi*10/48;

B2 = 2.5*conv([1,-2*ro1*cos(oo1),ro1^2],[1,-2*ro2*cos(oo2),ro2^2]);

%Poles
rp1 = 0.95;
op1 = 2*pi*2/48;

A2 = [1 -rp1*2*cos(op1) rp1^2];

%% Cascaded
BT = 25*conv(B1,B2);
AT = conv(A1,A2);

Output = filter(BT,AT,Input);

%Magnitude and Phase Spectra
F = [100:10:10000];
Fs = 48000;

figure
freqz(BT,AT,F,Fs);
ax = findall(gcf, 'Type', 'axes');
set(ax, 'XScale', 'log');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELEC3104: Project - Cochlear Signal Processing
% Authors:  Jared Agius     z5113699
%           Kioni Ndirangu  z5111826
% Date:     27/05/18

% Discrete Time Fourier Transform Function from Moodle-TutorialLab1 Solutions

function [X1] = dtft( y, nfft)

ylen = length(y); 
w=linspace(-pi,pi,nfft);
X1=0;

    for n=1:ylen    
         X1= X1+y(n).*exp(-1i*n*w); 
    end

end


