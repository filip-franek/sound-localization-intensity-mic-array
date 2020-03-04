% Multi_SI_Proc.m Multiple sound source localization using 3D Sound
% Intensity Method and Blind Source Separation
%
%     Sound sources location estimation based on 3D sound intensity 
%     measurement at presence of 4 sound sources simultaneously playing.
%     Sound recordings are taken from an anechoic chamber.
%     
%      variables to change
%          'varName' :	description
%                       can have such and such values
%                        * 'var1': description
%                        * 'var2': description
% Date: 03/03/2020
% Author: Filip Franek, AudibleBits
%
% Copyright (C) 2020 Filip Franek
%%
% Simulation:
%  - plane wave is assumed, distances between microphones are calculated
% Fs........ sampling frequency
% azwDeg/elwDeg ..... choose angle of N sources
% sepFold.... location of data saved(/to be saved)
clear all; clc; close all;
meaDate = '270514';

% Add script's folder to Matlab path
root = mfilename('fullpath');
[rootDir,name,ext] = fileparts(root);
cd(rootDir);

% cd(['', meaDate])
mixFile = 'sp1234mix_0dB';
roomExp = 'anCham';
tabName = [meaDate,'_', roomExp];
% azwDeg = [-12, 55, -35, 30];
% elwDeg = [16, -16, -49, -11]; 

azwDeg = [-12.2, 54.4, 29.2, -34.8 ]; % anCham
elwDeg = [11.8, -15.3, -10.7, -45.5]; 
% azwDeg = [-12, 55 30 -35]; % acLab
% elwDeg = [16, -16, -11 -49];
matFile=cell(4,1);
matFile{1} = 'sp1_speech1_n10db';
matFile{2} = 'sp2_norm_n10db';
matFile{3} = 'sp3_speech3_n10db';
matFile{4} = 'sp4_speech4_n10db';
bkFile = 'bkNoise';
d=50*0.001;         % set distance between microphones
Fs = 2^13;              % audio sampling freq.
saveDat = 0;        % save all data?
calOpt = 0;         % phase calibartion? if real data..
% txtMat = 'cali1234'; calFold = 'matTrun';
normOpt = 0;        % normalization of separated sources?
% titleOpt = 1;       % delete title for printing
% set(0,'DefaultFigureVisible','off')
% make speech always first for further comparison
sigTyp = {'Speech1', 'Normal', 'Speech3', 'Speech4'};   % 'Normal', 'Laplac', 'Uniform'
secL = 1.5;     % secL always larger then secLx
secLx = 1;
samLen = Fs*secL;
samLenx = Fs*secLx;
sigNum = length(sigTyp);         % number of sources
c=343;rho=1.21; R = d/4*sqrt(6);
df=2^(0);       % frequency resolution in Hz (2^n, n = (-3:1:+6))
bw=100;         % 100Hz bandwidth over which is averaged later for one version
nfft = Fs/df;   % number of FFT points
%% Plot options
% M..multi, S..sigle, X..mix, d..detail, 0..origianl, 
angMdW = 1.5;   angMdWo = 2;    angMdS = '-';      % W..lineWidth, S..lineStyle
angSW = .2;     angSWo = 2.5;   angSS = '-';      % W..lineWidth, S..lineStyle
angMW = 1.0;    angMWo = 0.5;     angMS = '-'; angMoS = '--';      % S..single,M..multiple,X..mixture,o..original angle
angXW = 1.5;    angXWo = 2;     angXS = '.';
azRan =52;              % angle range for multiple sources
elRan =52;              % angle range (overwritten (az,el) for solid range sometimes)
sgRan = 90;             % angle range for single unmixed source
sgRanS = 90;            % angle range for single separated sources : ylim +-X degrees
mixRan = 90;            % angle range of mixed not separated signals
xOver = 0;            % plotting x below and over ef.bw
efBW = [500, 1500];     % effective freq. bw. of intensity array
yAx = [efBW(1)-xOver, efBW(2)+xOver]; % xlimit for frequency
labFontSize = 10;
legFontSize = 8;
titFontSize = 10;
legFontSizeX = 10;  % because mixed seems different to others
tfW = 1;	% time-frequency plot, width
tW = 1;
angDistS = 'Normal';
angDistM = 'Normal';    % choose distribution evaluation for 
% angDistS = 'tLocationScale';
% angDistM = 'tLocationScale';    % choose distribution evaluation for 

sepFold = 'matSep';
colorLine = distinguishable_colors(sigNum*10);
% figure; image(reshape(colorLine,[1 size(colorLine)]))
%% 1:mix 2:single 3:time 4:multi 5:multi detail 6:single only 7: Error
plotAnMSTMSSo = [0 0 1 1 0 0 1]; % plot angle?
%%%%%%%%%%%%%%%% 1 2 3 4 5 6 7 8%%%%%%%%%%%%%%%%%%%
saveTab = 0;
saveEps = 0;        % save figures which are plotted
saveFig = 0;
prob = 0;           % plot probility distribution
simul = 0;          % simulation or measurement
genBkN = 1;         % add bk noise to simulated data
overWrite = 0;      % overwrite     !!!!!!!!?
plot3d_multi = 0;
spscOrg = [0 0 0];  % s...save, p...plot, sc.... soundsc original separated data, single data
% if wanna plot, plotAnMSTMSSo(6) == 0;
%% -----------SINGLE main-------------
%% LOAD Single Data(simul)
if simul == 0
matFold = '../Data/matOrg/';
for s = 1:sigNum    % load variables for s sources
    matPath = [matFold,matFile{s},'.mat'];
    load(matPath, 'php', 'Php')

    onset = 0;% in sec
    lenS = 1;% in sec

    phpTrun = php(onset*Fs+1:(onset+lenS)*Fs,:)';
    % figure; plot(phpTrun(1,:))
    % create single vectors of numMic pressure vectors
    for n = 1:size(php,2)
        phpTmp = phpTrun(n,:);
        eval(sprintf('p%dsg(s,:) = phpTmp;', n));
    end
	bkPath = [matFold, bkFile,'.mat'];

    load(bkPath);
    bk1s = bk(1:Fs);
    [BK, f, x, fh2] = anAudio(bk1s, Fs, 0);
    BKspl = abs(20*log10((BK)/2e-5));
end
    % load calibration data
    calFile1 = 'cali1234_norm_n55dB';
    matPath1 = [matFold,calFile1,'.mat'];
    load(matPath1);
    phpTrun = php(onset*Fs+1:(onset+lenS)*Fs,:)';
    for n = 1:4
        phpTmp = phpTrun(n,:);
        eval(sprintf('h1p%d = phpTmp;', n));
    end
    calFile2 = 'cali4567_norm_n55dB';
    matPath2 = [matFold,calFile2,'.mat'];
    load(matPath2);
    phpTrun = php(onset*Fs+1:(onset+lenS)*Fs,:)';
    for n = 4:7
        phpTmp = phpTrun(n,:);
        eval(sprintf('h2p%d = phpTmp;', n));
    end
elseif simul ==1
% adjusting the time delay for each microphone
% assuming plane wave p = A*e^(1*w*t)
% azwDeg = [  -12,  55, 30, -35];        % azimuth angle in deg <+180;-180>
% elwDeg = [   12, -15.3, -10.7, -45.5];        % elevation angle in deg <+180;-180>
angTxt = ['sim140614'];
sepName = ['SEP',angTxt,'.mat'];
sepData = [sepFold, '\', sepName];

r= 1000;                % distance centroid of probe - source [m] (should not really matter since assuming plane wave)
Fm = 2^21;      % frequency for manipulating signal, max 2^21


% if exist(sepData, 'file')&&overWrite==1
%     delete(sepData)
% end
% % check for number of signals
% if sigNum~=1 % if 1 source, no separation and jump out to 'simulPlaneWsingle.m'
% % coordinate position of each microphone
% if (exist(sepData,'file')==0) || (overWrite==1)
% simulation of sound source for convolutive mixtures
%                 X                Y                 Z
m1_pos = [ (-d)/(2*sqrt(3)),    (-d)/2,      (-d)/(2*sqrt(6)) ]; % mic1 at bottom back right
m2_pos = [ (-d)/(2*sqrt(3)),     (d)/2,      (-d)/(2*sqrt(6)) ]; % mic2 at bottom back left
m3_pos = [    (d)/(sqrt(3)),         0,      (-d)/(2*sqrt(6)) ]; % mic3 at bottom front
m4_pos = [                0,         0,     (d*(sqrt(6)))/(4) ]; % mic4 at top vertex
m5_pos = [  (d)/(2*sqrt(3)),    (-d)/2,      (-d)/(2*sqrt(6)) ]; % mic5 at bottom front right
m6_pos = [  (d)/(2*sqrt(3)),     (d)/2,      (-d)/(2*sqrt(6)) ]; % mic6 at bottom front left
m7_pos = [   -(d)/(sqrt(3)),         0,      (-d)/(2*sqrt(6)) ]; % mic7 at bottom back

M = [m1_pos; m2_pos; m3_pos; m4_pos; m5_pos; m6_pos; m7_pos];  % Matrix M missing first column of ones
x_m1 = M(1,1); y_m1 = M(1,2); z_m1 = M(1,3); m1_3d = [x_m1, y_m1, z_m1];
x_m2 = M(2,1); y_m2 = M(2,2); z_m2 = M(2,3); m2_3d = [x_m2, y_m2, z_m2];
x_m3 = M(3,1); y_m3 = M(3,2); z_m3 = M(3,3); m3_3d = [x_m3, y_m3, z_m3];
x_m4 = M(4,1); y_m4 = M(4,2); z_m4 = M(4,3); m4_3d = [x_m4, y_m4, z_m4]; 
x_m5 = M(5,1); y_m5 = M(5,2); z_m5 = M(5,3); m5_3d = [x_m5, y_m5, z_m5];
x_m6 = M(6,1); y_m6 = M(6,2); z_m6 = M(6,3); m6_3d = [x_m6, y_m6, z_m6];
x_m7 = M(7,1); y_m7 = M(7,2); z_m7 = M(7,3); m7_3d = [x_m7, y_m7, z_m7];

micNum = size(M,1);
p1res = zeros(1, secLx*Fm); p2res = zeros(1, secLx*Fm); p3res = zeros(1, secLx*Fm); p4res = zeros(1, secLx*Fm);
p5res = zeros(1, secLx*Fm); p6res = zeros(1, secLx*Fm); p7res = zeros(1, secLx*Fm);
T = samLen;      % #samples (Tmax)
t = (0:T-1)/Fs;         % time[s]
Tv = (0:T-1);           % sample vector
c = 343;                % speed of sound 20°
for s=1:sigNum      % simulate pressure propagation for each source..
    azwRad = deg2rad(azwDeg(s));
    elwRad = -deg2rad(elwDeg(s));   % because otherwise opposite sign for some reason..!!?!?

    % assuming spherical wave p = A/r*e^(1*w*t)
    z_s1 = sin(elwRad)*r;
    xy_s1 = cos(elwRad)*r;
    y_s1 = sin(azwRad)*xy_s1;
    x_s1 = cos(azwRad)*xy_s1;   % s1_3d = [x_s1; y_s1; z_s1];

    s1_m1 = sqrt((x_s1-x_m1)^2 + (y_s1-y_m1)^2 + (z_s1-z_m1)^2);
    s1_m2 = sqrt((x_s1-x_m2)^2 + (y_s1-y_m2)^2 + (z_s1-z_m2)^2);
    s1_m3 = sqrt((x_s1-x_m3)^2 + (y_s1-y_m3)^2 + (z_s1-z_m3)^2);
    s1_m4 = sqrt((x_s1-x_m4)^2 + (y_s1-y_m4)^2 + (z_s1-z_m4)^2);
    s1_m5 = sqrt((x_s1-x_m5)^2 + (y_s1-y_m5)^2 + (z_s1-z_m5)^2);
    s1_m6 = sqrt((x_s1-x_m6)^2 + (y_s1-y_m6)^2 + (z_s1-z_m6)^2);
    s1_m7 = sqrt((x_s1-x_m7)^2 + (y_s1-y_m7)^2 + (z_s1-z_m7)^2);
    
    Tm = 1/Fm;
    td = r/c;   % td_ind = round(td/Tm);
    Tsm1 = s1_m1/c;     Tsm1_ind = round(Tsm1/Tm);
    Tsm2 = s1_m2/c;     Tsm2_ind = round(Tsm2/Tm);
    Tsm3 = s1_m3/c;     Tsm3_ind = round(Tsm3/Tm);
    Tsm4 = s1_m4/c;     Tsm4_ind = round(Tsm4/Tm);    
    Tsm5 = s1_m5/c;     Tsm5_ind = round(Tsm5/Tm);
    Tsm6 = s1_m6/c;     Tsm6_ind = round(Tsm6/Tm);
    Tsm7 = s1_m7/c;     Tsm7_ind = round(Tsm7/Tm);
%     Tsm = [Tsm1_ind, Tsm2_ind, Tsm3_ind, Tsm4_ind, Tsm5_ind, Tsm6_ind, Tsm7_ind];
    minTsm = min([Tsm1_ind, Tsm2_ind, Tsm3_ind, Tsm4_ind, Tsm5_ind, Tsm6_ind, Tsm7_ind]);
    Tsm0(s,1:micNum) = [Tsm1_ind, Tsm2_ind, Tsm3_ind, Tsm4_ind, Tsm5_ind, Tsm6_ind, Tsm7_ind]-minTsm;

    % signal generator
    strSw = char(sigTyp(s));
    sigOrg = zeros(1, samLen);
    % scaling constant
    scNorm =  0.8368; scUni = 1.73; scLap = 0.7;
    scSp1 = 30; scSp2 = 15; scSp3 = 41; scSp4 = 12.2;
    
    switch strSw
        case 'Normal'   % normal distribution
        sigOrg= randn(samLen, 1)*scNorm; % 2 second - Gaussian distribution (amp <-5,5>)
        sigRes(s,:) = resample(sigOrg,Fm,Fs);  % resample for delay manipulation
        
        case 'Uniform'  % uniform distribution from -1 to 1
        sigOrg  = (1-2*rand(samLen, 1))*scUni; % 2 second amplitude (amp <-0.5,0.5>)
        sigRes(s,:) = resample(sigOrg,Fm,Fs);  % resample for delay manipulation
        
        case 'Laplac'   % Laplace distribution
        sigOrg  = random('exp',1,1,samLen)*scLap-random('exp',1,1,samLen)*scLap;
        sigRes(s,:) = resample(sigOrg,Fm,Fs);  % resample for delay manipulation
        
        case 'Speech1'  % female speech 1
        [voice Fsp] = audioread('wav\voiceOrg2s.wav');
        if size(voice)<(Fsp*secL)
            voicePad = padarray(voice,[(Fsp*secL-size(voice,1)), 0], 0, 'post');   % start array with p1, [a b]=> a=number of zeroes in column
            sigOrg(1:(Fsp*secL)) =  (scSp1*voicePad(1:(Fsp*secL)))';
        else
            sigOrg(1:(Fsp*secL)) =  (scSp1*voice(1:(Fsp*secL)))';
        end
        % original signal up-resample for delay manipulation
        sigTmp(s,:) = resample(sigOrg,Fm,Fsp);
        sigRes(s,:) = sigTmp(s,1:Fm*secL);
        
        case 'Speech2'  % male speech 2 (kennedy with bk noise)
        [voice Fsk] = audioread('wav\voice2Org2s.wav');
        if size(voice)<(Fsk*secL)
            voicePad = padarray(voice,[(Fsk*secL-size(voice,1)), 0], 0, 'post');   % start array with p1, [a b]=> a=number of zeroes in column
            sigOrg(1:(Fsk*secL)) =  (scSp2*voicePad(1:(Fsk*secL)))';
        else
            sigOrg(1:(Fsk*secL)) =  (scSp2*voice(1:(Fsk*secL)))';
        end
        % original signal up-resample for delay manipulation
        sigTmp(s,:) = resample(sigOrg,Fm,Fsk);
        sigRes(s,:) = sigTmp(s,1:Fm*secL);
        
        case 'Speech3'
        [voice Fsk] = audioread('wav\speech3Org2s.wav');
        if size(voice)<(Fsk*secL)
            voicePad = padarray(voice,[(Fsk*secL-size(voice,1)), 0], 0, 'post');   % start array with p1, [a b]=> a=number of zeroes in column
            sigOrg(1:(Fsk*secL)) =  (scSp3*voicePad(1:(Fsk*secL)))';
        else
            sigOrg(1:(Fsk*secL)) =  (scSp3*voice(1:(Fsk*secL)))';
        end
        % original signal up-resample for delay manipulation
        sigTmp(s,:) = resample(sigOrg,Fm,Fsk);
        sigRes(s,:) = sigTmp(s,1:Fm*secL);
        
        case 'Speech4'
        [voice Fsk] = audioread('wav\speech4org2s.wav');
        if size(voice)<(Fsk*secL)
            voicePad = padarray(voice,[(Fsk*secL-size(voice,1)), 0], 0, 'post');   % start array with p1, [a b]=> a=number of zeroes in column
            sigOrg(1:(Fsk*secL)) =  (scSp4*voicePad(1:(Fsk*secL)))';
        else
            sigOrg(1:(Fsk*secL)) =  (scSp4*voice(1:(Fsk*secL)))';
        end
        % original signal up-resample for delay manipulation
        sigTmp(s,:) = resample(sigOrg,Fm,Fsk);
        sigRes(s,:) = sigTmp(s,1:Fm*secL);
        
        otherwise
            disp('wrong input argument for sigGen')
    end % end of switch

    xG = sigRes(s,:);
    % generated signal time shifted
    xG1(s,:) = sigRes(s,Tsm0(s,1)+1:(secLx*Fm+Tsm0(s,1)));
    xG2(s,:) = sigRes(s,Tsm0(s,2)+1:(secLx*Fm+Tsm0(s,2)));
    xG3(s,:) = sigRes(s,Tsm0(s,3)+1:(secLx*Fm+Tsm0(s,3)));
    xG4(s,:) = sigRes(s,Tsm0(s,4)+1:(secLx*Fm+Tsm0(s,4)));
    xG5(s,:) = sigRes(s,Tsm0(s,5)+1:(secLx*Fm+Tsm0(s,5)));
    xG6(s,:) = sigRes(s,Tsm0(s,6)+1:(secLx*Fm+Tsm0(s,6)));
    xG7(s,:) = sigRes(s,Tsm0(s,7)+1:(secLx*Fm+Tsm0(s,7)));
    
    % generated delayed single signal down-sampled
    p1sgwoN(s,:) = resample(xG1(s,:),Fs,Fm);
    p2sgwoN(s,:) = resample(xG2(s,:),Fs,Fm);
    p3sgwoN(s,:) = resample(xG3(s,:),Fs,Fm);
    p4sgwoN(s,:) = resample(xG4(s,:),Fs,Fm);
    p5sgwoN(s,:) = resample(xG5(s,:),Fs,Fm);
    p6sgwoN(s,:) = resample(xG6(s,:),Fs,Fm);
    p7sgwoN(s,:) = resample(xG7(s,:),Fs,Fm);    
    % generated delayed signals mixed
    p1res = p1res + xG(Tsm0(s,1)+1:(secLx*Fm+Tsm0(s,1))); %/s1_m1; % for sph. wave
    p2res = p2res + xG(Tsm0(s,2)+1:(secLx*Fm+Tsm0(s,2))); %/s1_m2;
    p3res = p3res + xG(Tsm0(s,3)+1:(secLx*Fm+Tsm0(s,3))); %/s1_m3;
    p4res = p4res + xG(Tsm0(s,4)+1:(secLx*Fm+Tsm0(s,4))); %/s1_m4;    
    p5res = p5res + xG(Tsm0(s,5)+1:(secLx*Fm+Tsm0(s,5))); %/s1_m1; % for sph. wave
    p6res = p6res + xG(Tsm0(s,6)+1:(secLx*Fm+Tsm0(s,6))); %/s1_m2;
    p7res = p7res + xG(Tsm0(s,7)+1:(secLx*Fm+Tsm0(s,7))); %/s1_m3;
end % end of sound pressure simulation
% generated delayed, mixed signal down-sampled
p1woN = resample(p1res,Fs,Fm);
p2woN = resample(p2res,Fs,Fm);
p3woN = resample(p3res,Fs,Fm);
p4woN = resample(p4res,Fs,Fm);
p5woN = resample(p5res,Fs,Fm);
p6woN = resample(p6res,Fs,Fm);
p7woN = resample(p7res,Fs,Fm);

% generate bk noise
% scBkN = 0.0172; %
% save('bkN_sim', 'bkN', 'Bk')
load('bkN_sim', 'bkN', 'Bk')
if genBkN == 1  % add artificial bk noise to the signal
    for m = 1:micNum    % load variables for s sources
%         bkN(m,:) = randn(1, Fs)*scBkN;
%         [Bk(m,:), f, t, hAn1] = anAudio(bkN(m,:), Fs, 0);
        eval(sprintf('p%dorg = p%dwoN+bkN(%d,:);', m,m,m));
        for s = 1:sigNum    % load variables for s sources
            eval(sprintf('p%dsg(%d,:) = p%dsgwoN(%d,:)+bkN(%d,:);', m,s,m,s,m));
            eval(sprintf('x(%d,:) = p1sgwoN(%d,:);', s,s));
            [X, f, t, hAn1] = anAudio(x(s,:), Fs, 0);
SNRsg(s) = 20*log10((sqrt(mean(abs(X(efBW(1):efBW(2))).^2))/sqrt(mean(abs(Bk(m,efBW(1):efBW(2)).^2)))));
        end
    end% ['p',num2str(m),'woN']
    
    [P1, f, t, hAn1] = anAudio(p1woN, Fs, 0);
    SNRband = 20*log10((sqrt(mean(abs(P1(efBW(1):efBW(2))).^2))/sqrt(mean(abs(Bk(1,efBW(1):efBW(2)).^2)))));
elseif genBkN == 0
    for m = 1:micNum    % load variables for s sources
        eval(sprintf('p%dorg = p%dwoN;', m,m));
        eval(sprintf('p%dsg = p%dsgwoN;', m,m));
    end
end
bk1s = bkN(1,:);
BK=Bk(1,:);
BKspl = abs(20*log10((BK)/2e-5));
end
%% COMPUTE single sources angle
for s = 1:sigNum
%% Calibration? 
if calOpt == 1
    % tfestimate : Txy = Pyx(f) / Pxx(f), transform function estimate.
    % input real=>output Txy is a column vector of length nfft/2+1 for nfft evenand (nfft+1)/2 for nfft odd
    % [Txy,F]  = tfestimate(x,y,window,noverlap,nfft,fs,'twosided')
    [H41 fout] = tfestimate(h1p4,h1p1,nfft*3/4,50,nfft,Fs,'onesided');  % 'onesided' from 0 to the Nyquist frequency
    [H42 fout] = tfestimate(h1p4,h1p2,nfft*3/4,50,nfft,Fs,'onesided');
    [H43 fout] = tfestimate(h1p4,h1p3,nfft*3/4,50,nfft,Fs,'onesided');

    [H45 fout] = tfestimate(h2p4,h2p5,nfft*3/4,50,nfft,Fs,'onesided');
    [H46 fout] = tfestimate(h2p4,h2p6,nfft*3/4,50,nfft,Fs,'onesided');
    [H47 fout] = tfestimate(h2p4,h2p7,nfft*3/4,50,nfft,Fs,'onesided');
%     [H48 fout] = tfestimate(h3p4,h3p8,nfft*3/4,50,nfft,Fs,'onesided');
    
    pm14=-angle(H41);pm24=-angle(H42);pm34=-angle(H43);
    pm54=-angle(H45);pm64=-angle(H46);pm74=-angle(H47);%pm84=-angle(H48);
end
    %% Intensity
    % module 1
    [G21 fout] = cpsd(p2sg(s,:),p1sg(s,:),nfft,0,nfft,Fs,'onesided');
    [G31 fout] = cpsd(p3sg(s,:),p1sg(s,:),nfft,0,nfft,Fs,'onesided');
    [G41 fout] = cpsd(p4sg(s,:),p1sg(s,:),nfft,0,nfft,Fs,'onesided');
    [G32 fout] = cpsd(p3sg(s,:),p2sg(s,:),nfft,0,nfft,Fs,'onesided');
    [G42 fout] = cpsd(p4sg(s,:),p2sg(s,:),nfft,0,nfft,Fs,'onesided');
    [G43 fout] = cpsd(p4sg(s,:),p3sg(s,:),nfft,0,nfft,Fs,'onesided');
    % module 2
    [G45 fout] = cpsd(p4sg(s,:),p5sg(s,:),nfft,0,nfft,Fs,'onesided');
    [G46 fout] = cpsd(p4sg(s,:),p6sg(s,:),nfft,0,nfft,Fs,'onesided');
    [G47 fout] = cpsd(p4sg(s,:),p7sg(s,:),nfft,0,nfft,Fs,'onesided');
    [G65 fout] = cpsd(p6sg(s,:),p5sg(s,:),nfft,0,nfft,Fs,'onesided');
    [G75 fout] = cpsd(p7sg(s,:),p5sg(s,:),nfft,0,nfft,Fs,'onesided');
    [G76 fout] = cpsd(p7sg(s,:),p6sg(s,:),nfft,0,nfft,Fs,'onesided');

    if calOpt == 1  % Phase calibration
    g21=G21.*exp(1j*( pm24-pm14));   
    g31=G31.*exp(1j*( pm34-pm14));   g32=G32.*exp(1j*( pm34-pm24));
    g41=G41.*exp(1j*(-pm14));        g42=G42.*exp(1j*(-pm24));       g43=G43.*exp(1j*(-pm34));
    g65=G65.*exp(1j*(pm64-pm54));   g75=G75.*exp(1j*( pm74-pm54)); g45=G45.*exp(1j*(-pm54));
    g76=G76.*exp(1j*(pm74-pm64));   g46=G46.*exp(1j*(-pm64));      g47=G47.*exp(1j*(-pm74));
    % ..... g8...
    else
        g21 = G21; g31 = G31; g32 = G32; g41 = G41; g42 = G42; g43 = G43;
        g45 = G45; g46 = G46; g47 = G47; g65 = G65; g75 = G75; g76 = G76;
    end
    freq=fout;
    
    % Intensity of single sources
    % Module 1
    IXC1=(3*g31+3*g32+g41+g42-2*g43)./(4*sqrt(3)*rho*2*pi.*freq*d);
    IYC1=(2*g21+g31+g41-g32-g42)./(4*rho*2*pi*freq*d);
    IZC1=(g41+g42+g43)./(sqrt(6)*rho*2*pi*freq*d);
    % module 2
    IXC2=(-3*g75-3*g76-g45-g46+2*g47)./(4*sqrt(3)*rho*2*pi.*freq*d);
    IYC2=(2*g65+g75+g45-g76-g46)./(4*rho*2*pi.*freq*d);
    IZC2=(g45+g46+g47)./(sqrt(6)*rho*2*pi*freq*d);
    % basically IX=imag(IXC), imaginary../1/(rho*w*d)?
    IXs1(s,:)=imag(IXC1);
    IYs1(s,:)=imag(IYC1);
    IZs1(s,:)=imag(IZC1);
    IXs2(s,:)=imag(IXC2);
    IYs2(s,:)=imag(IYC2);
    IZs2(s,:)=imag(IZC2);
    IXs(s,:)=IXs1(s,:)+IXs2(s,:);
    IYs(s,:)=IYs1(s,:)+IYs2(s,:);
    IZs(s,:)=IZs1(s,:)+IZs2(s,:);
%% Angle
    if simul == 0
        azNs(s,:)=atan2(IYs(s,:),IXs(s,:))*(180/pi);

% for i=1:length(IXs(s,:)) % phase estimation probe 1
%     if IXs(s,i)<0 % azimuth for Ix<0
%         azNs(s,i)=atan2(IYs(s,i),IXs(s,i))*(180/pi)+180;
%     else        % azimuth for Ix>0
%         azNs(s,i)=atan2(IYs(s,i),IXs(s,i))*(180/pi);
%     end
% end
%  elNs(s,:)=atan2(IZs(s,:),sqrt(IYs(s,:).^2+IXs(s,:).^2)).*(180/pi); % elevation
    end
    if simul == 1;
        azNs(s,:)=atan2(IYs(s,:),IXs(s,:))*(180/pi)-180;
    end
    elNs(s,:)=(atan2(IZs(s,:),sqrt(IYs(s,:).^2+IXs(s,:).^2)).*(180/pi))';
    for n=1:length(azNs(s,:)) % + cluster similar angles to solid angle +-180 around truth angle
        if azNs(s,n) > azwDeg(s)+180
            azNsv(s,n) = azNs(s,n)-360;
        elseif azNs(s,n) < azwDeg(s)-180
            azNsv(s,n) = azNs(s,n)+360;
        else
            azNsv(s,n) = azNs(s,n);
        end   
        if elNs(s,n) > elwDeg(s)+180
            elNsv(s,n) = elNs(s,n)-360;
        elseif elNs(s,n) < elwDeg(s)-180
            elNsv(s,n) = elNs(s,n)+360;
        else
            elNsv(s,n) = elNs(s,n);
        end
    end
    
% statistical evaluation
    if simul == 0   % simulation returns angles which are shifted by pi/2
    azNs1(s,:)=atan2(IYs1(s,:),IXs1(s,:))*(180/pi);
    elseif simul ==1
	azNs1(s,:)=atan2(IYs1(s,:),IXs1(s,:))*(180/pi)-180;
    end
    elNs1(s,:)=(atan2(IZs1(s,:),sqrt(IYs1(s,:).^2+IXs1(s,:).^2)).*(180/pi))';
   
    for n=1:length(azNs1(s,:)) % + cluster similar angles to solid angle +-180 around truth angle
        if azNs1(s,n) > azwDeg(s)+180
            azNs1v(s,n) = azNs1(s,n)-360;
        elseif azNs1(s,n) < azwDeg(s)-180
            azNs1v(s,n) = azNs1(s,n)+360;
        else
            azNs1v(s,n) = azNs1(s,n);
        end   
        if elNs1(s,n) > elwDeg(s)+180
            elNs1v(s,n) = elNs1(s,n)-360;
        elseif elNs1(s,n) < elwDeg(s)-180
            elNs1v(s,n) = elNs1(s,n)+360;
        else
            elNs1v(s,n) = elNs1(s,n);
        end
    end
    if simul ==0
    azNs2(s,:)=atan2(IYs2(s,:),IXs2(s,:))*(180/pi);
    elseif simul ==1
    azNs2(s,:)=atan2(IYs2(s,:),IXs2(s,:))*(180/pi)-180;    
    end
    elNs2(s,:)=(atan2(IZs2(s,:),sqrt(IYs2(s,:).^2+IXs2(s,:).^2)).*(180/pi))';
    for n=1:length(azNs2(s,:)) % + cluster similar angles to solid angle +-180 around truth angle
        if azNs2(s,n) > azwDeg(s)+180
            azNs2v(s,n) = azNs2(s,n)-360;
        elseif azNs2(s,n) < azwDeg(s)-180
            azNs2v(s,n) = azNs2(s,n)+360;
        else
            azNs2v(s,n) = azNs2(s,n);
        end   
        if elNs2(s,n) > elwDeg(s)+180
            elNs2v(s,n) = elNs2(s,n)-360;
        elseif elNs2(s,n) < elwDeg(s)-180
            elNs2v(s,n) = elNs2(s,n)+360;
        else
            elNs2v(s,n) = elNs2(s,n);
        end
    end
    % azNmv(s,:) = azNm; elNmv(s,:) = elNm; % if pure calculated results
% 1st and 2nd order statistics
    azNs_arit(s) = roundn(sum(azNsv(s,efBW(1):efBW(2)))/length(azNsv(s,efBW(1):efBW(2))),-1);   % mean(azNsv(s,efBW(1):efBW(2)))
    elNs_arit(s) = roundn(sum(elNsv(s,efBW(1):efBW(2)))/length(elNsv(s,efBW(1):efBW(2))),-1);
    azNs1_arit(s) = roundn(sum(azNs1v(s,efBW(1):efBW(2)))/length(azNs1v(s,efBW(1):efBW(2))),-1);
    elNs1_arit(s) = roundn(sum(elNs1v(s,efBW(1):efBW(2)))/length(elNs1v(s,efBW(1):efBW(2))),-1);
    azNs2_arit(s) = roundn(sum(azNs2v(s,efBW(1):efBW(2)))/length(azNs2v(s,efBW(1):efBW(2))),-1);
    elNs2_arit(s) = roundn(sum(elNs2v(s,efBW(1):efBW(2)))/length(elNs2v(s,efBW(1):efBW(2))),-1);

    azNs_std(s) = roundn(sqrt(sum((azNsv(s,efBW(1):efBW(2))-azNs_arit(s)).^2)/length(azNsv(s,efBW(1):efBW(2)))), -2);    % std(azNmv(s,efBW(1):efBW(2)), 1)
    elNs_std(s) = roundn(sqrt(sum((elNsv(s,efBW(1):efBW(2))-elNs_arit(s)).^2)/length(elNsv(s,efBW(1):efBW(2)))), -2);
    azNs1_std(s) = roundn(sqrt(sum((azNs1v(s,efBW(1):efBW(2))-azNs1_arit(s)).^2)/length(azNs1v(s,efBW(1):efBW(2)))), -2);
    elNs1_std(s) = roundn(sqrt(sum((elNs1v(s,efBW(1):efBW(2))-elNs1_arit(s)).^2)/length(elNs1v(s,efBW(1):efBW(2)))), -2);
    azNs2_std(s) = roundn(sqrt(sum((azNs2v(s,efBW(1):efBW(2))-azNs2_arit(s)).^2)/length(azNs2v(s,efBW(1):efBW(2)))), -2);
    elNs2_std(s) = roundn(sqrt(sum((elNs2v(s,efBW(1):efBW(2))-elNs2_arit(s)).^2)/length(elNs2v(s,efBW(1):efBW(2)))), -2);
end
%% PLOT angle of single sources: detail
if plotAnMSTMSSo(2)==1
for s = 1:sigNum
    hFig70(s) = figure;  % azimuth angle (narrow band)
    hFs70(s) = subplot(1,2,1);
    h70(s) = plot(freq,azNsv(s,:),angSS, 'LineWidth',angSW,'Color',colorLine(s,:));
%      uistack(h70(s),'top')
    hold on
    h71(s) = plot(freq,azwDeg(s)*ones(length(freq), 1),'--','LineWidth',angSWo,'Color',colorLine(s,:)*.5);
    h72(s) = plot(freq,azNs1v(s,:),angSS, 'LineWidth',angSW,'Color',colorLine(19,:));
    h73(s) = plot(freq,azNs2v(s,:),angSS, 'LineWidth',angSW,'Color',colorLine(20,:));
    xlim(yAx);
    ylim([-sgRanS+azwDeg(s) sgRanS+azwDeg(s)]);
%     ylim([-90 90]);
    legAz1 = ['twi \theta_',num2str(s),'=',num2str( azNs_arit(s)),'°'];
    legAz2 = ['org \theta_',num2str(s),'=',num2str(azwDeg(s)),'°'];
    legAz3 = ['sg_1 \theta_',num2str(s),'=',num2str( azNs1_arit(s)),'°'];
    legendAz4 = ['sg_2 \theta_',num2str(s),'=',num2str( azNs2_arit(s)),'°'];    
    title_azN1 = ['Azimuth angle est. of ',num2str(s),'. source alone'];
    hT70(s) = title(title_azN1,'fontsize',titFontSize);
    xlabel('Frequency (Hz)','fontsize',labFontSize);
    ylabel('Angle (°)','fontsize',labFontSize);
    a = 'Ignoring extra legend entries.';
    MSGID ='MATLAB:legend:IgnoringExtraEntries';
    warning('off', MSGID)
    l70(s) = legend(legAz1,legAz2,legAz3,legendAz4, 'fontsize',legFontSize);
    uistack(h70(s),'top');
    uistack(h71(s),'top');
    set(l70(s),'Location','EastOutside');
    
    %         setgrid('xy','yx'); % create grid and shaded areas. grid on does not function. due to this suppresed fake warning
    x = [0 0 efBW(1) efBW(1) 0]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
    h1 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
    x = [efBW(2) efBW(2) efBW(2)+efBW(1) efBW(2)+efBW(1) efBW(2)]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
    h2 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');

    hFs80(s) = subplot(1,2,2);
    h80(s) = plot(freq,elNsv(s,:),angSS, 'LineWidth',angSW,'Color',colorLine(s,:));
    hold on
    h81(s) = plot(freq,elwDeg(s)*ones(length(freq), 1),'--','LineWidth',angSWo,'Color',colorLine(s,:)*.5);
    h82(s) = plot(freq,elNs1v(s,:),angSS, 'LineWidth',angSW,'Color',colorLine(s+4,:));
    h83(s) = plot(freq,elNs2v(s,:),angSS, 'LineWidth',angSW,'Color',colorLine(s+5,:));        
    xlim(yAx)
    ylim([-sgRanS+elwDeg(s) sgRanS+elwDeg(s)]);
%     ylim([-90 90]);
    legEl1 = ['twi \phi_',num2str(s),'=',num2str(elNs_arit(s)),'°'];
    legEl2 = ['org \phi_',num2str(s),'=',num2str(elwDeg(s)),'°'];
    legEl3 = ['sg_1 \phi',num2str(s),'=',num2str( elNs1_arit(s)),'°'];
    legendEl4 = ['sg_2 \phi',num2str(s),'=',num2str( elNs2_arit(s)),'°'];    
    title_azN1 = ['Elevation angle est. of ',num2str(s),'. source alone'];
    hT80(s) = title(title_azN1,'fontsize',titFontSize );
    xlabel('Frequency (Hz)','fontsize',labFontSize )
    ylabel('Angle (°)','fontsize',labFontSize )
    a = 'Ignoring extra legend entries.';
    MSGID ='MATLAB:legend:IgnoringExtraEntries';
    warning('off', MSGID)
    l80(s) = legend(legEl1,legEl2,legEl3,legendEl4, 'fontsize',legFontSize );
    uistack(h80(s),'top');
    uistack(h81(s),'top');
    set(l80(s),'Location','EastOutside');
    
    %         setgrid('xy','yx');
    x = [0 0 efBW(1) efBW(1) 0]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
    h1 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
    x = [efBW(2) efBW(2) efBW(2)+efBW(1) efBW(2)+efBW(1) efBW(2)]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
    h2 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
    set(hFig70(s), 'Position', [ 70 300 900 200]);
    ax1=get(hFs70(s),'Position'); % [0.1300 0.1900 0.0829 0.6950];
    ax1 = [0.1300 0.1900 0.0700 0.6950];
    ax2=get(hFs80(s),'Position'); % [0.5703 0.1900 0.0829 0.6950];
    ax2 = [0.5703 0.1900 0.0700 0.6950];
    % [shifts sideways, shifts up-down, stretches rightwards, stretches up]
    set(hFs70(s),'Position',[ax1(1)-.06, ax1(2), ax1(3)*3.5, ax1(4)]);
    set(hFs80(s),'Position',[ax2(1)-.02, ax2(2), ax2(3)*3.5, ax2(4)]);
    lx1=get(l70(s),'Position'); % [0.1300, 0.5838, 0.6512, 0.3412]
    set(l70(s) ,'Position', [lx1(1)-0.01, lx1(2), lx1(3), lx1(4)]);%
    lx1=get(l80(s),'Position'); % [0.1300, 0.5838, 0.6512, 0.3412]
    set(l80(s) ,'Position', [lx1(1)-0.01, lx1(2), lx1(3), lx1(4)]);%
end
end
%% ++++++++++++MULTI main+++++++++++++
%% skip if only single
if plotAnMSTMSSo(6)==0
%% LOAD Multi Data
if simul == 0;
matFold = '../Data/matOrg/';
matPath = [matFold,mixFile,'.mat'];
load(matPath)

onset = 0;% in sec
lenS = 1;% in sec
phpTrun = php(onset*Fs+1:(onset+lenS)*Fs,:)';
% figure; plot(phpTrun(1,:))    
% create single vectors of pressure vectors
for n = 1:size(tmp,2)
%     pTmp = tmpTrun(n,:);
%     eval(sprintf('p%dorg = pTmp;', n));
    phpTmp = phpTrun(n,:);  % save filtered signal into pXorg
    eval(sprintf('p%dorg = phpTmp;', n));
end
end
% forming microphone group
x1234 = [p1org;p2org;p3org;p4org];  x1234Name = varname(x1234);
x4567 = [p4org;p5org;p6org;p7org];  x4567Name = varname(x4567);
x1234567 = [p1org;p2org;p3org;p4org;p5org;p6org;p7org];  x1234567Name = varname(x1234567);
% x12345678 = [p1org;p2org;p3org;p4org;p5org;p6org;p7org;p8org];  x12345678Name = varname(x12345678);
if overWrite == 1
% for ICA, projection doesn't work, GCC-PHAT the best with Toeplitz and
% hierarchical clustering
% try hierarchical, pitch shifted, but maybe better. toeplitz is good
% BARBI, GCC_PHAT, Fuzzy c-means seems better (don't choose normal)
% BARBI, GCC_PHAT, hierarchical, binary seems OK as well
% 
if simul == 1
    ICAmet = 'barbi';     % bgl efica extefica bwasobi
    onSet = 0;
    simil = 2;      %1...projections 2...GCC-PHAT 3...Coherence
    clust = 'rfcm';     % cluster: rfcm=fuzzy c-mean, hclus=hierarchical
    weiTyp = 5;         % varies-1:normal,2:toeplitz,3: binary...
    fiLen = 20;         % filter length
%     ICAmet = 'bgl';     % bgl efica extefica bwasobi
%     onSet = 0;
%     simil = 1;      %1...projections 2...GCC-PHAT 3...Coherence
%     clust = 'rfcm';     % cluster: rfcm=fuzzy c-mean, hclus=hierarchical
%     weiTyp = 1;         % varies-1:normal,2:toeplitz,3: binary...
%     fiLen = 20;         % filter length
elseif simul == 0
    ICAmet = 'efica';     % bgl efica extefica bwasobi
    onSet = 0;
    simil = 'gcc';      %1...projections 2...GCC-PHAT 3...Coherence
    subBandsR = 1;
    muPar = 1;
    clust = 'hclus';     % cluster: rfcm=fuzzy c-mean, hclus=hierarchical
    weiTyp = 'toep';         % varies-1:normal,2:toeplitz,3: binary...
    weiVal = 1;    
    fiLen = 20;         % filter length
end
[y1, est1, data1, Xdis1] = BASS_tabcd(x1234(:,onSet+1:Fs+onSet),Fs, ICAmet, fiLen, simil, subBandsR, muPar, clust, weiTyp, weiVal, x1234Name);
[y2, est2, data2, Xdis2] = BASS_tabcd(x4567(:,onSet+1:Fs+onSet),Fs, ICAmet, fiLen, simil, subBandsR, muPar, clust, weiTyp, weiVal, x4567Name);
% [y1, est1, data1] = BBS_tabcd(x1234(:,onSet+1:Fs+onSet),Fs, ICAmet, 1, simil, clust, weiTyp, fiLen, x1234Name);
% [y2, est2, data2] = BBS_tabcd(x4567(:,onSet+1:Fs+onSet),Fs, ICAmet, 1, simil, clust, weiTyp, fiLen, x4567Name);
% [y, est, data] = BBS_tabcd(x12345678(:,onSet+1:Fs+onSet),Fs, ICAmet, 1, simil, clust, weiTyp, fiLen, x12345678Name);
% [y, est, data] = BASS_tabcd(x1234567(:,onSet+1:Fs+onSet),Fs, ICAmet, fiLen, simil, subBandsR, muPar, clust, weiTyp, weiVal, x1234567Name);
% [y, est, data] = BASS_tabcd(x1234567,Fs, 'bgl', 20, gcc, 1, 1, 'rfcm', 'norm', 1, x1234567Name);
onPost = 0; % post processing onset due to peaks at the start

p1 = zeros(sigNum, size(est1,2));
p2 = zeros(sigNum, size(est1,2));
p3 = zeros(sigNum, size(est1,2));
p4 = zeros(sigNum, size(est1,2));
p4b = zeros(sigNum, size(est1,2));
p5 = zeros(sigNum, size(est1,2));
p6 = zeros(sigNum, size(est1,2));
p7 = zeros(sigNum, size(est1,2));
% p8 = zeros(sigNum, size(est1,2));
%% compare and sort out signals
% for k = 1:sigNum
%     % est_resp(i,:,j) ith microphone, jth source
%     p1Tmp(k,:) = est1(1,onPost+1:Fs+onPost,k);
%     p2Tmp(k,:) = est1(2,onPost+1:Fs+onPost,k);
%     p3Tmp(k,:) = est1(3,onPost+1:Fs+onPost,k);
%     p4Tmp(k,:) = est1(4,onPost+1:Fs+onPost,k);
%     
%     p4bTmp(k,:)= est2(1,onPost+1:Fs+onPost,k);
%     p5Tmp(k,:) = est2(2,onPost+1:Fs+onPost,k);
%     p6Tmp(k,:) = est2(3,onPost+1:Fs+onPost,k);
%     p7Tmp(k,:) = est2(4,onPost+1:Fs+onPost,k);
% % 	p8(k,:) = est1(4,onPost+1:Fs+onPost,k);   
% % similarity measure between 
%     c1(:,k) = mscohere(p1sg(1,:),p1Tmp(k,:));  
%     c1val(k) = sum(c1(:,k))/length(c1(:,k));
%     c2(:,k) = mscohere(p1sg(2,:),p1Tmp(k,:));
%     c2val(k) = sum(c2(:,k))/length(c2(:,k));
%     c3(:,k) = mscohere(p1sg(3,:),p1Tmp(k,:));
%     c3val(k) = sum(c3(:,k))/length(c3(:,k));
%     c4(:,k) = mscohere(p1sg(4,:),p1Tmp(k,:));
%     c4val(k) = sum(c4(:,k))/length(c4(:,k));
%     
%     c4b(:,k) = mscohere(p1sg(1,:),p4bTmp(k,:));
%     c4bval(k) = sum(c4b(:,k))/length(c4b(:,k));
%     c5(:,k) = mscohere(p1sg(2,:),p4bTmp(k,:));  % similarity measure between 
%     c5val(k) = sum(c5(:,k))/length(c5(:,k));
%     c6(:,k) = mscohere(p1sg(3,:),p4bTmp(k,:));
%     c6val(k) = sum(c6(:,k))/length(c6(:,k));
%     c7(:,k) = mscohere(p1sg(4,:),p4bTmp(k,:));
%     c7val(k) = sum(c7(:,k))/length(c7(:,k));
% end
% % swap channels as they were originally for correcting angle evaluation
% [c1max ind1]= max(c1val);
% [c2max ind2]= max(c2val);
% [c3max ind3]= max(c3val);
% [c4max ind4]= max(c4val);
% 
% [c4bmax ind4b]=max(c4bval);
% [c5max ind5]= max(c5val);
% [c6max ind6]= max(c6val);
% [c7max ind7]= max(c7val);
% 
% indC1 = [ind1, ind2, ind3, ind4];
% indC2 = [ind4b, ind5, ind6, ind7];
% 
% for k = 1:sigNum
%     p1(k,:) = p1Tmp(indC1(k),:);
%     p2(k,:) = p2Tmp(indC1(k),:);
%     p3(k,:) = p3Tmp(indC1(k),:);
%     p4(k,:) = p4Tmp(indC1(k),:);
%     
%     p4b(k,:)= p4bTmp(indC2(k),:);
%     p5(k,:) = p5Tmp(indC2(k),:);
%     p6(k,:) = p6Tmp(indC2(k),:);
%     p7(k,:) = p7Tmp(indC2(k),:);    
% end

% manual manipulation
% indC1 = [2, 3, 4, 1];
% indC2 = [1, 3, 4, 2];
% % FOR     'bwasobi','pro','rfcm','norm'
% indC1 = [1, 4, 2, 3];% [1, 3, 4, 2]
% indC2 = [1, 4, 2, 3];
% FOR     'efica','pro','rfcm','norm'
indC1 = [2, 4, 3, 1]; % [4, 1, 3, 2]
indC2 = [3, 4, 1, 2]; % [3, 4, 1, 2]
% FOR     'efica','gcc','hclus','norm'
indC1 = [4, 2, 3, 1]; % [4, 2, 3, 1]
indC2 = [1, 2, 4, 3]; % [1, 2, 4, 3]

% indC1 = [4, 3, 2, 1];   % for by seeing
% indC2 = [3, 4, 1, 2];
% indC1 = [3 4 1 2];   % for simulated data bgl
% indC2 = [3 4 1 2];   % 2, 1, 4, 3
% indC1 = [1 2 3 4];   % for simulated data efica
% indC2 = [2 1 4 3];   % 2, 1, 4, 3
for k = 1:sigNum    % kth signal
    % est_resp(i,:,j) ith microphone, jth source
    p1Tmp(k,:) = est1(1,onPost+1:Fs+onPost,k);
    p2Tmp(k,:) = est1(2,onPost+1:Fs+onPost,k);
    p3Tmp(k,:) = est1(3,onPost+1:Fs+onPost,k);
    p4Tmp(k,:) = est1(4,onPost+1:Fs+onPost,k);
    
    p4bTmp(k,:) = est2(1,onPost+1:Fs+onPost,k);
    p5Tmp(k,:) = est2(2,onPost+1:Fs+onPost,k);
    p6Tmp(k,:) = est2(3,onPost+1:Fs+onPost,k);
    p7Tmp(k,:) = est2(4,onPost+1:Fs+onPost,k);    
end
for k = 1:sigNum
    p1(k,:) = p1Tmp(indC1(k),:);
    p2(k,:) = p2Tmp(indC1(k),:);
    p3(k,:) = p3Tmp(indC1(k),:);
    p4(k,:) = p4Tmp(indC1(k),:);
    
    p4b(k,:) = p4bTmp(indC2(k),:);
    p5(k,:) = p5Tmp(indC2(k),:);
    p6(k,:) = p6Tmp(indC2(k),:);
    p7(k,:) = p7Tmp(indC2(k),:);    
end
elseif overWrite == 0;
if simul == 0
load('../Data/matSep\SEPanec270514.mat', 'p1', 'p2', 'p3', 'p4', 'p4b', 'p5', 'p6', 'p7')
    %     load('matSep\SEP030614.mat', 'p1', 'p2', 'p3', 'p4', 'p4b', 'p5', 'p6', 'p7')
    %     load('matSep\SEPlab150614.mat', 'p1', 'p2', 'p3', 'p4', 'p4b', 'p5', 'p6', 'p7')
%     save('matSep\SEPanec270514.mat', 'p1', 'p2', 'p3', 'p4', 'p4b', 'p5', 'p6', 'p7')
elseif simul ==1
%     load('matSep\SEPsim120614.mat', 'p1', 'p2', 'p3', 'p4', 'p4b', 'p5', 'p6', 'p7')
    load('../Data/matSep\SEPsim150614.mat', 'p1', 'p2', 'p3', 'p4', 'p4b', 'p5', 'p6', 'p7')
%     save('matSep\SEPsim120614.mat', 'p1', 'p2', 'p3', 'p4', 'p4b', 'p5', 'p6', 'p7')
end
end
%% COMPUTE multi source intensity
% CPSD (Cross power spectral density) ~ used to be CSD
% Calculates Single-Sided Cross Spectral Density Gnm(w)
% [Pxy,F] = cpsd(x,y,window,noverlap,nfft,fs, '...sided')
% estimates the cross powerspectral density Pxy of the discrete-time signals x and y usingthe Welch's averaged, modified periodogram method of spectral estimation.The cross power spectral density is the distribution of power perunit frequency and is defined as
% For real x and y, cpsd returns a one-sided CPSD and for complex x or y,it returns a two-sided CPSD.
azMax = max(azwDeg); azMin = min(azwDeg); elMax = max(elwDeg); elMin = min(elwDeg); 
for s = 1:sigNum
% Intensity (multiple)
%     % normalize data
%     p1n(s,:) = p1(s,:)/max(abs(p1(s,:)));
%     p2n(s,:) = p2(s,:)/max(abs(p2(s,:)));
%     p3n(s,:) = p3(s,:)/max(abs(p3(s,:)));
%     p4n(s,:) = p1(s,:)/max(abs(p1(s,:)));

    if normOpt ==1;
    % module 1
    [G21 fout] = cpsd(p2n(s,:),p1n(s,:),nfft,0,nfft,Fs,'onesided');
    [G31 fout] = cpsd(p3n(s,:),p1n(s,:),nfft,0,nfft,Fs,'onesided');
    [G41 fout] = cpsd(p4n(s,:),p1n(s,:),nfft,0,nfft,Fs,'onesided');
    [G32 fout] = cpsd(p3n(s,:),p2n(s,:),nfft,0,nfft,Fs,'onesided');
    [G42 fout] = cpsd(p4n(s,:),p2n(s,:),nfft,0,nfft,Fs,'onesided');
    [G43 fout] = cpsd(p4n(s,:),p3n(s,:),nfft,0,nfft,Fs,'onesided');
    % module 2
    [G45 fout] = cpsd(p4n(s,:),p5n(s,:),nfft,0,nfft,Fs,'onesided');
    [G46 fout] = cpsd(p4n(s,:),p6n(s,:),nfft,0,nfft,Fs,'onesided');
    [G47 fout] = cpsd(p4n(s,:),p7n(s,:),nfft,0,nfft,Fs,'onesided');
    [G65 fout] = cpsd(p6n(s,:),p5n(s,:),nfft,0,nfft,Fs,'onesided');
    [G75 fout] = cpsd(p7n(s,:),p5n(s,:),nfft,0,nfft,Fs,'onesided');
    [G76 fout] = cpsd(p7n(s,:),p6n(s,:),nfft,0,nfft,Fs,'onesided');    
    else
    % module 1
    [G21 fout] = cpsd(p2(s,:),p1(s,:),nfft,0,nfft,Fs,'onesided');
    [G31 fout] = cpsd(p3(s,:),p1(s,:),nfft,0,nfft,Fs,'onesided');
    [G41 fout] = cpsd(p4(s,:),p1(s,:),nfft,0,nfft,Fs,'onesided');
    [G32 fout] = cpsd(p3(s,:),p2(s,:),nfft,0,nfft,Fs,'onesided');
    [G42 fout] = cpsd(p4(s,:),p2(s,:),nfft,0,nfft,Fs,'onesided');
    [G43 fout] = cpsd(p4(s,:),p3(s,:),nfft,0,nfft,Fs,'onesided');
    % module 2
    [G45 fout] = cpsd(p4b(s,:),p5(s,:),nfft,0,nfft,Fs,'onesided');
    [G46 fout] = cpsd(p4b(s,:),p6(s,:),nfft,0,nfft,Fs,'onesided');
    [G47 fout] = cpsd(p4b(s,:),p7(s,:),nfft,0,nfft,Fs,'onesided');
    [G65 fout] = cpsd(p6(s,:),p5(s,:),nfft,0,nfft,Fs,'onesided');
    [G75 fout] = cpsd(p7(s,:),p5(s,:),nfft,0,nfft,Fs,'onesided');
    [G76 fout] = cpsd(p7(s,:),p6(s,:),nfft,0,nfft,Fs,'onesided');    
    end

    if calOpt == 1  % Phase calibration
    % tfestimate : Txy = Pyx(f) / Pxx(f), transform function estimate.
    % input real=>output Txy is a column vector of length nfft/2+1 for nfft evenand (nfft+1)/2 for nfft odd
    % [Txy,F]  = tfestimate(x,y,window,noverlap,nfft,fs,'twosided')

    pm14=-angle(H41);pm24=-angle(H42);pm34=-angle(H43);
    pm54=-angle(H45);pm64=-angle(H46);pm74=-angle(H47);% pm84=-angle(H48);
    g21=G21.*exp(1j*( pm24-pm14));   
    g31=G31.*exp(1j*( pm34-pm14));   g32=G32.*exp(1j*( pm34-pm24));
    g41=G41.*exp(1j*(-pm14));        g42=G42.*exp(1j*(-pm24));       g43=G43.*exp(1j*(-pm34));
    g65=G65.*exp(1j*(pm64-pm54));   g75=G75.*exp(1j*( pm74-pm54)); g45=G45.*exp(1j*(-pm54));
    g76=G76.*exp(1j*(pm74-pm64));   g46=G46.*exp(1j*(-pm64));      g47=G47.*exp(1j*(-pm74));
    % ..... g8...
    else
        g21 = G21; g31 = G31; g32 = G32; g41 = G41; g42 = G42; g43 = G43;
        g45 = G45; g46 = G46; g47 = G47; g65 = G65; g75 = G75; g76 = G76;
    end
    freq=fout;
    % Complex intensity
    % Module 1
    IXC1=(3*g31+3*g32+g41+g42-2*g43)./(4*sqrt(3)*rho*2*pi.*freq*d);
    IYC1=(2*g21+g31+g41-g32-g42)./(4*rho*2*pi*freq*d);
    IZC1=(g41+g42+g43)./(sqrt(6)*rho*2*pi*freq*d);
    % module 2
    IXC2=(-3*g75-3*g76-g45-g46+2*g47)./(4*sqrt(3)*rho*2*pi.*freq*d);
    IYC2=(2*g65+g75+g45-g76-g46)./(4*rho*2*pi.*freq*d);
    IZC2=(g45+g46+g47)./(sqrt(6)*rho*2*pi*freq*d);
    % basically IX=imag(IXC), imaginary../1/(rho*w*d)?
    IXm1(s,:)=imag(IXC1);
    IYm1(s,:)=imag(IYC1);
    IZm1(s,:)=imag(IZC1);
    IXm2(s,:)=imag(IXC2);
    IYm2(s,:)=imag(IYC2);
    IZm2(s,:)=imag(IZC2);
    IXm(s,:)=IXm1(s,:)+IXm2(s,:);
    IYm(s,:)=IYm1(s,:)+IYm2(s,:);
    IZm(s,:)=IZm1(s,:)+IZm2(s,:);
%% COMPUTE multi source angle
% Double twisted tetrahedron
    if simul == 1;
        azNm(s,:)=atan2(IYm(s,:),IXm(s,:))*(180/pi)-180;
    elseif simul ==0
        azNm(s,:)=atan2(IYm(s,:),IXm(s,:))*(180/pi);
    end
    elNm(s,:)=(atan2(IZm(s,:),sqrt(IYm(s,:).^2+IXm(s,:).^2)).*(180/pi))';
    for n=1:length(azNm(s,:)) % + cluster similar angles to solid angle +-180 around truth angle
        if azNm(s,n) > azwDeg(s)+180
            azNmv(s,n) = azNm(s,n)-360;
        elseif azNm(s,n) < azwDeg(s)-180
            azNmv(s,n) = azNm(s,n)+360;
        else
            azNmv(s,n) = azNm(s,n);
        end   
        if elNm(s,n) > elwDeg(s)+180
            elNmv(s,n) = elNm(s,n)-360;
        elseif elNm(s,n) < elwDeg(s)-180
            elNmv(s,n) = elNm(s,n)+360;
        else
            elNmv(s,n) = elNm(s,n);
        end
    end
% Single tetrahedron 1
    if simul == 1;
    azNm1(s,:)=atan2(IYm1(s,:),IXm1(s,:))*(180/pi)+180;
    elseif simul == 0;
        azNm1(s,:)=atan2(IYm1(s,:),IXm1(s,:))*(180/pi);
    end
    elNm1(s,:)=(atan2(IZm1(s,:),sqrt(IYm1(s,:).^2+IXm1(s,:).^2)).*(180/pi))';
    for n=1:length(azNm1(s,:)) % + cluster similar angles to solid angle +-180 around truth angle
        if azNm1(s,n) > azwDeg(s)+180
            azNm1v(s,n) = azNm1(s,n)-360;
        elseif azNm1(s,n) < azwDeg(s)-180
            azNm1v(s,n) = azNm1(s,n)+360;
        else
            azNm1v(s,n) = azNm1(s,n);
        end   
        if elNm1(s,n) > elwDeg(s)+180
            elNm1v(s,n) = elNm1(s,n)-360;
        elseif elNm1(s,n) < elwDeg(s)-180
            elNm1v(s,n) = elNm1(s,n)+360;
        else
            elNm1v(s,n) = elNm1(s,n);
        end
    end
% Single tetrahedron 2
    if simul == 1
    azNm2(s,:)=atan2(IYm2(s,:),IXm2(s,:))*(180/pi)+180;
    elseif simul == 0
        azNm2(s,:)=atan2(IYm2(s,:),IXm2(s,:))*(180/pi);
    end
    elNm2(s,:)=(atan2(IZm2(s,:),sqrt(IYm2(s,:).^2+IXm2(s,:).^2)).*(180/pi))';
    for n=1:length(azNm2(s,:)) % + cluster similar angles to solid angle +-180 around truth angle
        if azNm2(s,n) > azwDeg(s)+180
            azNm2v(s,n) = azNm2(s,n)-360;
        elseif azNm2(s,n) < azwDeg(s)-180
            azNm2v(s,n) = azNm2(s,n)+360;
        else
            azNm2v(s,n) = azNm2(s,n);
        end   
        if elNm2(s,n) > elwDeg(s)+180
            elNm2v(s,n) = elNm2(s,n)-360;
        elseif elNm2(s,n) < elwDeg(s)-180
            elNm2v(s,n) = elNm2(s,n)+360;
        else
            elNm2v(s,n) = elNm2(s,n);
        end
    end
    % azNmv(s,:) = azNm; elNmv(s,:) = elNm; % if pure calculated results

% 1st and 2nd order statistics
    azNm_arit(s) = roundn(sum(azNmv(s,efBW(1):efBW(2)))/length(azNmv(s,efBW(1):efBW(2))),-1);   % mean(azNmv(s,efBW(1):efBW(2)))
    elNm_arit(s) = roundn(sum(elNmv(s,efBW(1):efBW(2)))/length(elNmv(s,efBW(1):efBW(2))),-1);
    azNm1_arit(s) = roundn(sum(azNm1v(s,efBW(1):efBW(2)))/length(azNm1v(s,efBW(1):efBW(2))),-1);
    elNm1_arit(s) = roundn(sum(elNm1v(s,efBW(1):efBW(2)))/length(elNm1v(s,efBW(1):efBW(2))),-1);
    azNm2_arit(s) = roundn(sum(azNm2v(s,efBW(1):efBW(2)))/length(azNm2v(s,efBW(1):efBW(2))),-1);
    elNm2_arit(s) = roundn(sum(elNm2v(s,efBW(1):efBW(2)))/length(elNm2v(s,efBW(1):efBW(2))),-1);

    azNm_std(s) = roundn(sqrt(sum((azNmv(s,efBW(1):efBW(2))-azNm_arit(s)).^2)/length(azNmv(s,efBW(1):efBW(2)))), -2);    % std(azNmv(s,efBW(1):efBW(2)), 1)
    elNm_std(s) = roundn(sqrt(sum((elNmv(s,efBW(1):efBW(2))-elNm_arit(s)).^2)/length(elNmv(s,efBW(1):efBW(2)))), -2);
    azNm1_std(s) = roundn(sqrt(sum((azNm1v(s,efBW(1):efBW(2))-azNm1_arit(s)).^2)/length(azNm1v(s,efBW(1):efBW(2)))), -2);
    elNm1_std(s) = roundn(sqrt(sum((elNm1v(s,efBW(1):efBW(2))-elNm1_arit(s)).^2)/length(elNm1v(s,efBW(1):efBW(2)))), -2);
    azNm2_std(s) = roundn(sqrt(sum((azNm2v(s,efBW(1):efBW(2))-azNm2_arit(s)).^2)/length(azNm2v(s,efBW(1):efBW(2)))), -2);
    elNm2_std(s) = roundn(sqrt(sum((elNm2v(s,efBW(1):efBW(2))-elNm2_arit(s)).^2)/length(elNm2v(s,efBW(1):efBW(2)))), -2);
    
    azNm_med(s) = roundn(median(azNmv(s,efBW(1):efBW(2))),-1);
    elNm_med(s) = roundn(median(elNmv(s,efBW(1):efBW(2))),-1);
end
%% PLOT angle of multiple sources
if plotAnMSTMSSo(4) == 1
    legAz1=cell(sigNum,1);   % initialize legend entries
    legAz2=cell(sigNum,1);
    legEl1=cell(sigNum,1);
    legEl2=cell(sigNum,1);
    angMWo = 1;
    angMW = 0.5;
    df = 1;
    dfMark = 50;
    close all
    markCell = {'o', '*', 's', '^'};
for s = 1:sigNum
    hFig10 = figure(100);  % azimuth angle (narrow band)
    subplot(2,1,1)
    h10 = plot(freq(1:df:end),smooth(azNmv(s,1:df:end)),angMS, 'LineWidth',angMW,'Color',[0 0 0]);
    hold on
    plot(freq(1:dfMark:end),smooth(azNmv(s,1:dfMark:end)),'.', 'LineWidth',angMW,'Color',[0 0 0], 'marker', markCell{s});
    h11 = plot(freq,azwDeg(s)*ones(length(freq),1),angMoS,'LineWidth',angMWo,'Color',[0.4,0.4,0.4]);
    xlim(yAx);
%     ylim([-azRan+azMin azRan+azMax]);
    ylim([-90 90]);
    legAz1{s} = ['twi \theta_',num2str(s),'=',num2str(azNm_arit(s)),'°'];
    legAz2{s} = ['org \theta_',num2str(s),'=',num2str(azwDeg(s)),'°'];
    title_azN1 = ['Azimuth angle est. of ',num2str(sigNum),' sources'];
%     hT10 = title(title_azN1,'fontsize',titFontSize );
    xlabel('Frequency (Hz)','fontsize',labFontSize );
    ylabel('Azimuth angle (°)','fontsize',labFontSize );
    set(hFig10, 'Position', [ 70 10 800 620]);
    if s ==sigNum   % legend for 2 sources
        a = 'Ignoring extra legend entries.';
        MSGID ='MATLAB:legend:IgnoringExtraEntries';
        warning('off', MSGID)
        l1 = legend(legAz1{1},legAz2{1},legAz1{2},legAz2{2},legAz1{3},legAz2{3},legAz1{4},legAz2{4}, 'fontsize',legFontSize );
        
        set(l1,'Location','EastOutside');
        %         setgrid('xy','yx'); % create grid and shaded areas. grid on does not function. due to this suppresed fake warning
        x = [0 0 efBW(1) efBW(1) 0]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
        h1 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
        x = [efBW(2) efBW(2) efBW(2)+efBW(1) efBW(2)+efBW(1) efBW(2)]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
        h2 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
    end
    legend off
%     hFig20 = figure(20);  % elevation angle (narrow band)
    subplot(2,1,2)
    h20 = plot(freq,smooth(elNmv(s,:)),angMS, 'LineWidth',angMW,'Color',[0 0 0]);
    hold on
	plot(freq(1:dfMark:end),smooth(elNmv(s,1:dfMark:end)),'.', 'LineWidth',angMW,'Color',[0 0 0], 'marker', markCell{s});
    h21 = plot(freq,elwDeg(s)*ones(length(freq), 1),angMoS,'LineWidth',angMWo,'Color',[0.4,0.4,0.4]);
    xlim(yAx)
%     ylim([-azRan+elMin azRan+elMax])
    ylim([-90 90]);
    legEl1{s} = ['twi \phi_',num2str(s),'=',num2str(elNm_arit(s)),'°'];
    legEl2{s} = ['org \phi_',num2str(s),'=',num2str(elwDeg(s)),'°'];
    title_azN1 = ['Elevation angle est. of ',num2str(s),' sources'];
%     hT20 =title(title_azN1,'fontsize',titFontSize );
    xlabel('Frequency (Hz)','fontsize',labFontSize )
    ylabel('Elevation angle (°)','fontsize',labFontSize )
    if s ==sigNum   % legend for 2 sources
        l2 = legend(legEl1{1},legEl2{1},legEl1{2},legEl2{2},legEl1{3},legEl2{3},legEl1{4},legEl2{4}, 'fontsize',legFontSize);
        set(l2,'Location','EastOutside');
        
        %         setgrid('xy','yx');
        x = [0 0 efBW(1) efBW(1) 0]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
        h1 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
        x = [efBW(2) efBW(2) efBW(2)+efBW(1) efBW(2)+efBW(1) efBW(2)]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
        h2 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
    end
    legend off
    angLim = 180;
    x3Val = (-angLim):.5:(angLim);  % angles across the probability density is calculacted
    y3Val = (-angLim):.5:(angLim);
    pd3Az(s,:)= fitdist(azNmv(s, efBW(1):efBW(2))','normal'); % estimate the angle by probability distribution
    pd3El(s,:)= fitdist(elNmv(s, efBW(1):efBW(2))','normal');
    pdf3Az(s,:) = pdf(pd3Az(s,:),x3Val);     % probability density function
    pdf3El(s,:) = pdf(pd3El(s,:),y3Val);
end
end
%% PLOT angle of separated sources: detail
if plotAnMSTMSSo(5) == 1
for s = 1:sigNum
    hFig50(s) = figure;  % azimuth angle (narrow band)
    hFs50(s) = subplot(1,2,1);
    h50(s) = plot(freq,azNmv(s,:),angMdS, 'LineWidth',angMdW,'Color',colorLine(s,:));
    hold on
    h51(s) = plot(freq,azwDeg(s)*ones(length(freq), 1),'--','LineWidth',angMdWo,'Color',colorLine(s,:)*.5);
    h52(s) = plot(freq,azNm1v(s,:),angMdS, 'LineWidth',angMdW,'Color',colorLine(s+1,:));
    h53(s) = plot(freq,azNm2v(s,:),angMdS, 'LineWidth',angMdW,'Color',colorLine(s+2,:));
    xlim(yAx);
    ylim([-sgRan+azwDeg(s) sgRan+azwDeg(s)]);
%     ylim([-90 90]);
    % needs to be changed for single or more sources
    legAz1 = ['twi \theta_',num2str(s),'=',num2str(azNm_arit(s)),'°'];
    legAz2 = ['org \theta_',num2str(s),'=',num2str(azwDeg(s)),'°'];
    legAz3 = ['sg_1 \theta_',num2str(s),'=',num2str(azNm1_arit(s)),'°'];
    legendAz4 = ['sg_2 \theta_',num2str(s),'=',num2str(azNm2_arit(s)),'°'];
    title_azN1 = ['Azimuth angle est. of ',num2str(s),'. separated source'];
    hT50(s) = title(title_azN1,'fontsize',titFontSize );
    xlabel('Frequency (Hz)','fontsize',labFontSize );
    ylabel('Angle (°)','fontsize',labFontSize );
    l50(s) = legend(legAz1,legAz2,legAz3,legendAz4, 'fontsize',legFontSize );
    uistack(h50(s),'top');
    uistack(h51(s),'top');
    set(l50(s),'Location','EastOutside');
    
    %         setgrid('xy','yx'); % create grid and shaded areas. grid on does not function. due to this suppresed fake warning
	x = [0 0 efBW(1) efBW(1) 0]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
    h1 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
    x = [efBW(2) efBW(2) efBW(2)+efBW(1) efBW(2)+efBW(1) efBW(2)]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
    h2 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');

    %     hFig20 = figure(20);  % elevation angle (narrow band)
    hFs60(s) = subplot(1,2,2);
    h60(s) = plot(freq,elNmv(s,:),angMdS, 'LineWidth',angMdW,'Color',colorLine(s,:));
    hold on
    h61(s) = plot(freq,elwDeg(s)*ones(length(freq), 1),'--','LineWidth',angXWo,'Color',colorLine(s,:)*.5);
    xlim(yAx)
    h62(s) = plot(freq,elNm1v(s,:),angMdS, 'LineWidth',angMdW,'Color',colorLine(s+1,:));
    h63(s) = plot(freq,elNm2v(s,:),angMdS, 'LineWidth',angMdW,'Color',colorLine(s+2,:));
    ylim([-sgRan+elwDeg(s) sgRan+elwDeg(s)]);
    legEl1 = ['twi \phi_',num2str(s),'=',num2str(elNm_arit(s)),'°'];
    legEl2 = ['org \phi_',num2str(s),'=',num2str(elwDeg(s)),'°'];
    legEl3 = ['sg_1 \phi_',num2str(s),'=',num2str(elNm1_arit(s)),'°'];
    legendEl4 = ['sg_2 \phi_',num2str(s),'=',num2str(elNm2_arit(s)),'°'];
    title_azN1 = ['Elevation angle est. of ',num2str(s),'. separated source'];
    hT60(s) = title(title_azN1,'fontsize',titFontSize );
    xlabel('Frequency (Hz)','fontsize',labFontSize )
    ylabel('Angle (°)','fontsize',labFontSize )
    l60(s) = legend(legEl1,legEl2,legEl3, legendEl4, 'fontsize',legFontSize );
    uistack(h60(s),'top');
    uistack(h61(s),'top');
    set(l60(s),'Location','EastOutside');
    
    %         setgrid('xy','yx');
    x = [0 0 efBW(1) efBW(1) 0]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
    h1 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
    x = [efBW(2) efBW(2) efBW(2)+efBW(1) efBW(2)+efBW(1) efBW(2)]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
    h2 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
    set(hFig50(s), 'Position', [ 70 300 900 200]);
    ax1=get(hFs50,'Position'); % [0.1300 0.1900 0.0829 0.6950];
    ax1 = [0.1300 0.1900 0.0700 0.6950];
    ax2=get(hFs60,'Position'); % [0.5703 0.1900 0.0829 0.6950];
    ax2 = [0.5703 0.1900 0.0700 0.6950];
    % [shifts sideways, shifts up-down, stretches rightwards, stretches up]
    set(hFs50(s),'Position',[ax1(1)-.06, ax1(2), ax1(3)*3.5, ax1(4)]);
    set(hFs60(s),'Position',[ax2(1)-.02, ax2(2), ax2(3)*3.5, ax2(4)]);
    lx1=get(l50(s),'Position'); % [0.1300, 0.5838, 0.6512, 0.3412]
    set(l50(s) ,'Position', [lx1(1)-0.01, lx1(2), lx1(3), lx1(4)]);%
    lx1=get(l60(s),'Position'); % [0.1300, 0.5838, 0.6512, 0.3412]
    set(l60(s) ,'Position', [lx1(1)-0.01, lx1(2), lx1(3), lx1(4)]);%
end
end
%% COMPUTE mixed multiple sources angle
% calibration if chosen
    [G21 fout] = cpsd(p2org(1:Fs),p1org(1:Fs),nfft,0,nfft,Fs,'onesided');
    [G31 fout] = cpsd(p3org(1:Fs),p1org(1:Fs),nfft,0,nfft,Fs,'onesided');
    [G41 fout] = cpsd(p4org(1:Fs),p1org(1:Fs),nfft,0,nfft,Fs,'onesided');
    [G32 fout] = cpsd(p3org(1:Fs),p2org(1:Fs),nfft,0,nfft,Fs,'onesided');
    [G42 fout] = cpsd(p4org(1:Fs),p2org(1:Fs),nfft,0,nfft,Fs,'onesided');
    [G43 fout] = cpsd(p4org(1:Fs),p3org(1:Fs),nfft,0,nfft,Fs,'onesided');
    
    if calOpt == 1  % Phase calibration
    % tfestimate : Txy = Pyx(f) / Pxx(f), transform function estimate.
    % input real=>output Txy is a column vector of length nfft/2+1 for nfft evenand (nfft+1)/2 for nfft odd
    % [Txy,F]  = tfestimate(x,y,window,noverlap,nfft,fs,'twosided')
    calFile = ['',calFold,'\',txtMat, '.mat'];
    load(calFile, 'p1fil', 'p2fil', 'p3fil', 'p4fil')
    h1p4 = p4fil; h1p1 = p1fil; h1p2 = p2fil; h1p3 = p3fil;
    nfft = 2^13;
    [H41 fout] = tfestimate(h1p4,h1p1,nfft*3/4,50,nfft,Fs,'onesided');  % 'onesided' from 0 to the Nyquist frequency
    [H42 fout] = tfestimate(h1p4,h1p2,nfft*3/4,50,nfft,Fs,'onesided');
    [H43 fout] = tfestimate(h1p4,h1p3,nfft*3/4,50,nfft,Fs,'onesided');

    % pm14=angle(H41);pm24=angle(H42);pm34=angle(H43);
    pm14=-angle(H41);pm24=-angle(H42);pm34=-angle(H43);
    g21=G21.*exp(1j*( pm24-pm14));   
    g31=G31.*exp(1j*( pm34-pm14));   g32=G32.*exp(1j*( pm34-pm24));
    g41=G41.*exp(1j*(-pm14));        g42=G42.*exp(1j*(-pm24));       g43=G43.*exp(1j*(-pm34));
    else
        g21 = G21; g31 = G31; g32 = G32; g41 = G41; g42 = G42; g43 = G43;
    end
    freq=fout; 
    
    % Intensitly
    % module 1
    IXC=(3*g31+3*g32+g41+g42-2*g43)./(4*sqrt(3)*rho*2*pi.*freq*d);
    IYC=(2*g21+g31+g41-g32-g42)./(4*rho*2*pi*freq*d);
    IZC=(g41+g42+g43)./(sqrt(6)*rho*2*pi*freq*d);
    % basically IX=imag(IXC), imaginary../1/(rho*w*d)?
    IXx=imag(IXC);
    IYx=imag(IYC);
    IZx=imag(IZC);

    if simul == 1
        azNx=atan2(IYx,IXx)*(180/pi)+180;
    elseif simul == 0
        azNx=atan2(IYx,IXx)*(180/pi);
    end
    elNx=(atan2(IZx,sqrt(IYx.^2+IXx.^2)).*(180/pi))';
    for n=1:length(azNx) % cluster similar angles to solid angle +-180 around truth angle
        if azNx(n) > azwDeg(s)+180
            azNxv(n) = azNx(n)-360;
        elseif azNx(n) < azwDeg(s)-180
            azNxv(n) = azNx(n)+360;
        else
            azNxv(n) = azNx(n);
        end   
        if elNx(n) > elwDeg(s)+180
            elNxv(n) = elNx(n)-360;
        elseif elNx(n) < elwDeg(s)-180
            elNxv(n) = elNx(n)+360;
        else
            elNxv(n) = elNx(n);
        end
    end
    % statistical evaluation
    azNx_arit = roundn(sum(azNxv(efBW(1):efBW(2)))/length(azNxv(efBW(1):efBW(2))),-1);
    azNx_med = roundn(median(azNxv(efBW(1):efBW(2))),-1);
    elNx_arit = roundn(sum(elNxv(efBW(1):efBW(2)))/length(elNxv(efBW(1):efBW(2))),-1);
    elNx_med = roundn(median(elNxv(efBW(1):efBW(2))),-1);
%% PLOT angle of mixed multiple sources
if plotAnMSTMSSo(1) == 1
hFig30 = figure;  % azimuth angle (narrow band)
hFs30 = subplot(1,2,1);
h30 = plot(freq,azNxv,angXS, 'LineWidth',angXW,'Color',colorLine(12,:)); hold on;
for s=1:sigNum
h31(s) = plot(freq,azwDeg(s)*ones(length(freq), 1),'--','LineWidth',angXWo,'Color',colorLine(s,:));
legAz2x{s} = ['org \theta_{',num2str(s),'}=',num2str(azwDeg(s)),'°'];
end
xlim(yAx);
ylim([-mixRan+azNx_arit azNx_arit+mixRan]);
legAz1x = ['est \theta =',num2str(azNx_arit),'°'];
title_azN1 = ['Azimuth angle est. of ',num2str(sigNum),' mixed sources'];
hT30 = title(title_azN1,'fontsize',titFontSize );
xlabel('Frequency (Hz)','fontsize',labFontSize );
ylabel('Angle (°)','fontsize',labFontSize );
l30 = legend([legAz1x,legAz2x], 'fontsize',legFontSizeX );
set(l30,'Location','EastOutside');

%         setgrid('xy','yx'); % create grid and shaded areas. grid on does not function. due to this suppresed fake warning
x = [0 0 efBW(1) efBW(1) 0]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
h1 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
x = [efBW(2) efBW(2) efBW(2)+efBW(1) efBW(2)+efBW(1) efBW(2)]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
h2 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');

%     hFig40 = figure;  % elevation angle (narrow band)
hFs40 = subplot(1,2,2);
h40 = plot(freq,elNxv,angXS, 'LineWidth',angXW,'Color',colorLine(12,:)); hold on;
for s=1:sigNum
h41(s) = plot(freq,elwDeg(s)*ones(length(freq), 1),'--','LineWidth',angMWo,'Color',colorLine(s,:));
legEl2x{s} = ['org \phi_{',num2str(s),'}=',num2str(elwDeg(s)),'°'];
end
xlim(yAx)
ylim([-mixRan+azNx_arit azNx_arit+mixRan]);
legEl1x = ['est \phi =',num2str(elNx_arit),'°'];
title_azN1 = ['Elevation angle est. of ',num2str(sigNum),' mixed sources'];
hT40 = title(title_azN1,'fontsize',titFontSize );
xlabel('Frequency (Hz)','fontsize',labFontSize )
ylabel('Angle (°)','fontsize',labFontSize )
l40 = legend([legEl1x,legEl2x], 'fontsize',legFontSizeX);
set(l40,'Location','EastOutside');

%         setgrid('xy','yx');
x = [0 0 efBW(1) efBW(1) 0]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
h1 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
x = [efBW(2) efBW(2) efBW(2)+efBW(1) efBW(2)+efBW(1) efBW(2)]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
h2 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
set(hFig30, 'Position', [ 70 10 900 200]);

ax1=get(hFs30,'Position'); % [0.1300 0.1900 0.0829 0.6950];
ax1 = [0.1300 0.1900 0.0700 0.6950];
ax2=get(hFs40,'Position'); % [0.5703 0.1900 0.0829 0.6950];
ax2 = [0.5703 0.1900 0.0700 0.6950];
% [shifts sideways, shifts up-down, stretches rightwards, stretches up]
set(hFs30,'Position',[ax1(1)-.07, ax1(2), ax1(3)*4, ax1(4)]);
set(hFs40,'Position',[ax2(1)-.02, ax2(2), ax2(3)*4, ax2(4)]);
lx1=get(l30,'Position'); % [0.1300, 0.5838, 0.6512, 0.3412]
set(l30 ,'Position', [lx1(1)-0.01, lx1(2), lx1(3), lx1(4)]);%
lx1=get(l40,'Position'); % [0.1300, 0.5838, 0.6512, 0.3412]
set(l40 ,'Position', [lx1(1)-0.01, lx1(2), lx1(3), lx1(4)]);%
end
%% Error evaluation
if plotAnMSTMSSo(7) == 1
%% plot histogram
hFig60 = figure(60);    % histogram multiple sources
histRng = 15;           % angle display range +-histRng
binRng = 0.1;           % angle resolution
hFig65 = figure(65);  % 3D histogram + 2D plot azimuth
hFig66 = figure(66);  % 3D histogram + 2D plot elevation

for s = 1:sigNum;
azHist(s,:) = (azwDeg(s)-180):binRng:(azwDeg(s)+180);       % set histogram vector
elHist(s,:) = (elwDeg(s)-180):binRng:(elwDeg(s)+180);
binCount(s,:) = histc(azNm(s,efBW(1):efBW(2)),azHist(s,:)); % compute histogram for azi, no plot
binCount(s+4,:) = histc(elNm(s,efBW(1):efBW(2)),elHist(s,:));   % compute histogram for ele, no plot
end
[histVal histInd] = max(max(binCount));     % find maximum value of histogram
for s = 1:sigNum;
%     x3Val = (-180):.5:(180);  % angles across the probability density is calculacted
%     y3Val = (-180):.5:(180);
%     pd3Az(s,:)= fitdist(azNsv(s, efBW(1):efBW(2))',angDistM); % estimate the angle by probability distribution
%     pd3El(s,:)= fitdist(elNsv(s, efBW(1):efBW(2))',angDistM);
%     pdf3Az(s,:) = pdf(pd3Az(s,:),x3Val);     % probability density function
%     pdf3El(s,:) = pdf(pd3El(s,:),y3Val);
%     [f,xi] = ksdensity(x)
    
% plot histogram of single and multiple source angles
    figure(hFig60)
    hSub61(s) = subplot(2,4,s);
    hB61(s) = bar(azHist(s,:), binCount(s,:), 'FaceColor', [0.4 0.4 0.4],'EdgeColor',[0.4 0.4 0.4]);
    hold on
    [azR azC] = find(round(azHist(s,:),3)==round(azwDeg(s)),3); azC = azC(1);
    azBar = zeros(1,length(azHist(s,:)));
    azBar(azC) = histVal/10;
    hBb61(s) = bar(azHist(s,:), azBar,'FaceColor', 'k','EdgeColor','k', 'LineWidth',1);
    hold off
%     hist(hSub61(s), azNm(s,efBW(1):efBW(2)),azHist(s,:));
    leg61 = ['twi \theta_',num2str(s),'=',num2str(azNm_arit(s)),'°'];
    xlabel(['Azimuth \theta_', num2str(s), '(°)'],'fontsize',labFontSize );
    ylabel('Frequency (bins)','fontsize',labFontSize );
%     l61(s) = legend(leg61, 'fontsize',legFontSize );
    xlim([azwDeg(s)-histRng azwDeg(s)+histRng])
    ylim(hSub61(s), [0 histVal ]);
%     set(hBb61(s),'BaseValue',-1); 
    
    hSub62(s) = subplot(2,4,s+4);
    hB62(s) = bar(elHist(s,:), binCount(s+4,:), 'FaceColor', [0.4 0.4 0.4],'EdgeColor',[0.4 0.4 0.4]);
    hold on
    [elR elC] = find(round(elHist(s,:),1)==round(elwDeg(s),1)); elC = elC(1);
    elBar = zeros(1,length(elHist(s,:)));
    elBar(elC) = histVal/10;
    hBb62(s) = bar(elHist(s,:), elBar, 'FaceColor', 'k','EdgeColor','k', 'LineWidth',1);
    hold off
%     hist(hSub62(s), elNm(s,efBW(1):efBW(2)),elHist);
    leg62 = ['org \phi_',num2str(s),'=',num2str(elNm_arit(s)),'°'];
    xlabel(['Elevation \theta_', num2str(s), '(°)'],'fontsize',labFontSize );
    ylabel('Frequency (bins)','fontsize',labFontSize );
%     l62(s) = legend(leg62, 'fontsize',legFontSize );
    xlim([elwDeg(s)-histRng elwDeg(s)+histRng])
    ylim(hSub62(s), [0 histVal ]);
%     title_azN1 = ['Azimuth angle est. of ',num2str(sigNum),' sources'];
%     hT10 = title(title_azN1,'fontsize',titFontSize );
    
    figure(hFig65)
%     hSub66(s) = subplot(2,4,s);
%     figure(hFig69)
%     hSub69(s) = subplot(1,4,s);
    binRng3d = 1;	% width for x and y bins
    data = [azNm(s,efBW(1):efBW(2));elNm(s,efBW(1):efBW(2))]';
    x = [azwDeg(s)-histRng:binRng3d:azwDeg(s)+histRng ];	% use same range for y bins, can be different
    y = [elwDeg(s)-histRng:binRng3d:elwDeg(s)+histRng ]; % range for x bins

%     hist3(data, 'edges', {x,y});
%     %     axis([elwDeg(s)-histRng elwDeg(s)+histRng  azwDeg(s)-histRng azwDeg(s)+histRng 0 histVal])
%     xlim(hSub66(s),[azwDeg(s)-histRng azwDeg(s)+histRng])
%     ylim(hSub66(s),[elwDeg(s)-histRng elwDeg(s)+histRng])
%     zlim(hSub66(s),[0 histVal ]);    
% 
%     xlabel(['\theta_', num2str(s), '(°)'],'fontsize',labFontSize );
%     ylabel(['\phi_', num2str(s), '(°)'],'fontsize',labFontSize );
%     zlabel('Frequency (bins)','fontsize',labFontSize );
% 
%     set(gcf, 'renderer', 'opengl');
%     set(get(gca,'child'), 'FaceColor', 'interp', 'EdgeColor','k', 'CDataMode', 'auto');
%     view
    
    hSub67(s) = subplot(2,2,s);
    n = hist3(data, 'edges', {x,y}); % Extract histogram
    colormap(flipud(gray))
    hHex(s) = pcolor(x,y,n');
    xlabel(['Azimuth \theta_', num2str(s), '(°)'],'fontsize',labFontSize );
    ylabel(['Elevation \phi_', num2str(s), '(°)'],'fontsize',labFontSize );
    axis square
        set(gcf, 'renderer', 'opengl');
    set(get(gca,'child'), 'FaceColor', 'interp', 'EdgeColor','w', 'CDataMode', 'auto');

   
%     colormap
end
% % set(l61,'Location','EastOutside');
%     set(h,'FaceColor',[0 .5 .5],'EdgeColor','w')
set(hFig65, 'Position', [ 70 325 900 420]);
set(hFig60, 'Position', [ 70 10 900 420]);
%%
if plot3d_multi == 1
hFig69 = figure(69);
angPlot = 15;
for s = 1:sigNum;
    %     angPdf = pdfAz(s,:)+pdfEl(s,:);
    angPdf(s,:,:) = pdf3El(s,:)'*pdf3Az(s,:);
    probMax = max(max(max((angPdf(:,:,:)))));
%     surf(x3Val,y3Val,angPdf)
    hS1prob3D(s) = subplot(2,4,s);
    mesh(x3Val,y3Val,squeeze(angPdf(s,:,:))), view(0,90)
    xlim(hS1prob3D(s),[azwDeg(s)-angPlot azwDeg(s)+angPlot])
    ylim(hS1prob3D(s),[elwDeg(s)-angPlot elwDeg(s)+angPlot])
%     axis([-angPlot angPlot -angPlot angPlot 0 0.01])
    axis square
    hold on
    hS2prob3D(s) = subplot(2,4,s+4);
     surf(x3Val,y3Val,squeeze(angPdf(s,:,:)));
	xlim(hS2prob3D(s),[azwDeg(s)-angPlot azwDeg(s)+angPlot])
    ylim(hS2prob3D(s),[elwDeg(s)-angPlot elwDeg(s)+angPlot])
    zlim(hS2prob3D(s),[0 probMax]); 
% 	xlim([azwDeg(s)-histRng azwDeg(s)+histRng])
%     ylim([elwDeg(s)-histRng elwDeg(s)+histRng])
%     zlim([0 histVal ]); 
    axis square
    shading flat
    hold on
end
set(hFig69, 'Position', [ 70 10 900 420]);
%%
hFig68 = figure(68);
angPlot = 70;
for s = 1:sigNum;
    %     angPdf = pdfAz(s,:)+pdfEl(s,:);
    angPdf(s,:,:) = pdf3El(s,:)'*pdf3Az(s,:);
    probMax = max(max(max((angPdf(:,:,:)))));
%     surf(x3Val,y3Val,angPdf)
    hS3prob3D(s) = subplot(1,2,1);
    mesh(x3Val,y3Val,squeeze(angPdf(s,:,:))), view(0,90)
    xlim(hS3prob3D(s),[-angPlot +angPlot])
    ylim(hS3prob3D(s),[-angPlot +angPlot])
    zlim(hS3prob3D(s),[0 histVal ]); 
%     axis([-angPlot angPlot -angPlot angPlot 0 0.01])
    axis square
    hold on
    hS2prob3D(s) = subplot(1,2,2);
    surf(x3Val,y3Val,squeeze(angPdf(s,:,:)));
% 	xlim(hS2prob3D(s),[-angPlot +angPlot])
%     ylim(hS2prob3D(s),[-angPlot +angPlot])
    zlim(hS2prob3D(s),[0 probMax]); 
% 	xlim([azwDeg(s)-histRng azwDeg(s)+histRng])
%     ylim([elwDeg(s)-histRng elwDeg(s)+histRng])
%     zlim([0 histVal ]); 
    axis square
    shading flat
    hold on
end
set(hFig68, 'Position', [ 70 10 900 420]);

end
end
% Angle error calculation
azNs_err = azNs_arit-azwDeg;
azNs1_err = azNs1_arit-azwDeg;
azNs2_err = azNs2_arit-azwDeg;
azNm_err = azNm_arit-azwDeg;
azNm1_err = azNm1_arit-azwDeg;
azNm2_err = azNm2_arit-azwDeg;
elNs_err = elNs_arit-elwDeg;
elNs1_err = elNs1_arit-elwDeg;
elNs2_err = elNs2_arit-elwDeg;
elNm_err = elNm_arit-elwDeg;
elNm1_err = elNm1_arit-elwDeg;
elNm2_err = elNm2_arit-elwDeg;
% save('matAn\angleEvalSGlab',...     % save statistical evaluation of single sources
% 'azNs_arit','elNs_arit','azNs1_arit','elNs1_arit','azNs2_arit','elNs2_arit',...
% 'azNs_std','elNs_std','azNs1_std','elNs1_std','azNs2_std','elNs2_std',...
% 'azNs_err', 'azNs1_err', 'azNs2_err', 'azNm_err', 'azNm1_err', 'azNm2_err', ...
% 'elNs_err', 'elNs1_err', 'elNs2_err', 'elNm_err', 'elNm1_err', 'elNm2_err')

for s = 1:sigNum
% SNR
    [X(s,:), f, t, hAn1] = anAudio(p1sg(s,:), Fs, 0);
    Xspl(s,:) = abs(20*log10((X(s,:))/2e-5));
    SNRband(s) = 20*log10((sqrt(mean(abs(X(s,efBW(1):efBW(2))).^2))/sqrt(mean(abs(BK(efBW(1):efBW(2)).^2)))));
% Multi-angle probability
	if prob == 1
    pdAz(s,:)= fitdist(azNmv(s, efBW(1):efBW(2))',angDistM); % estimate the angle by probability distribution
    pdEl(s,:)= fitdist(elNmv(s, efBW(1):efBW(2))',angDistM); % estimate the angle by probability distribution
    xVal(s,:) = (-20+(azwDeg(s))):0.01:(20+(azwDeg(s)));  % angles across the probability density is calculacted
    yVal(s,:) = (-20+(elwDeg(s))):0.01:(20+(elwDeg(s)));  % angles across the probability density is calculacted
    pdfAz(s,:) = pdf(pdAz(s,:),xVal(s,:));     % probability density function
    pdfEl(s,:) = pdf(pdEl(s,:),yVal(s,:));     % probability density function
    dotL = 1;

    hFig6(s) = figure;
    h6a = plot(xVal(s,:),pdfAz(s,:),'b','LineWidth',1.5); hold on
    h6c = plot(yVal(s,:),pdfEl(s,:),'r', 'LineWidth',1.5);
    h6b = plot(azwDeg(s), dotL, 'k');
    h6d = plot(elwDeg(s), dotL, 'k');
    [max_p_az, index_p_az] = max(pdfAz(s,:));
    [max_p_el, index_p_el] = max(pdfEl(s,:));
    u_az(s) = xVal(s,index_p_az);
    u_el(s) = yVal(s,index_p_el);
    axis([min([xVal(s,:),yVal(s,:)]), max([xVal(s,:),yVal(s,:)]), 0 max([max_p_az,max_p_el])*1.1])
    
    legAz1 = ['mode \theta_',num2str(s),'=',num2str(u_az(s)),'°'];
    legAz2 = ['org \theta_',num2str(s),'=',num2str(azwDeg(s)),'°'];
    legEl1 = ['mode \phi_',num2str(s),'=',num2str(u_el(s)),'°'];
    legEl2 = ['org \phi_',num2str(s),'=',num2str(elwDeg(s)),'°'];
    l6 = legend(legAz1, legEl1, legAz2, legEl2);
    set(l6, 'Location', 'EastOutside')
    title1 = ['Angle est. of ',num2str(s),'. source'];
    title(title1,'fontsize',titFontSize );
    xlabel('Angle (°)','fontsize',labFontSize );
    ylabel('Probability (-)','fontsize',labFontSize );
    set(hFig6, 'Position', [ 70 10 700 180]);

    hFig7(s) = figure;
    hist(azNmv(s, efBW(1):efBW(2)),3600)
    title('distribution of ...')
    xlabel('magnitude [-]'); ylabel('[frequency]')
    end
% single source probability

    pdAz(s,:)= fitdist(azNsv(s, efBW(1):efBW(2))',angDistS); % estimate the angle by probability distribution
    pdEl(s,:)= fitdist(elNsv(s, efBW(1):efBW(2))',angDistS);
    xVal(s,:) = (-180+(azwDeg(s))):0.1:(180+(azwDeg(s)));  % angles across the probability density is calculacted
    yVal(s,:) = (-180+(elwDeg(s))):0.1:(180+(elwDeg(s)));
    pdfAz(s,:) = pdf(pdAz(s,:),xVal(s,:));     % probability density function
    pdfEl(s,:) = pdf(pdEl(s,:),yVal(s,:));
    
    if prob == 1
    hFig8(s) = figure;
    h8a = plot(xVal(s,:),pdfAz(s,:),'b','LineWidth',1.5); hold on
    h8c = plot(yVal(s,:),pdfEl(s,:),'r', 'LineWidth',1.5);
    %# vertical line
    h8b = graph2d.constantline(azwDeg(s), 'LineStyle',':', 'Color',[.7 .7 .7]);
    changedependvar(h8b,'x');
    h8d = graph2d.constantline(elwDeg(s), 'LineStyle',':', 'Color',[.5 .5 .5]);
    changedependvar(h8d,'x');
    [max_p_az, index_p_az] = max(pdfAz(s,:));
    [max_p_el, index_p_el] = max(pdfEl(s,:));
    u_az(s) = xVal(s,index_p_az);
    u_el(s) = yVal(s,index_p_el);
    axis([min([xVal(s,:),yVal(s,:)]), max([xVal(s,:),yVal(s,:)]), 0 max([max_p_az,max_p_el])*1.1])
    
    legAz1 = ['mode \theta_',num2str(s),'=',num2str(u_az(s)),'°'];
    legAz2 = ['org \theta_',num2str(s),'=',num2str(azwDeg(s)),'°'];
    legEl1 = ['mode \phi_',num2str(s),'=',num2str(u_el(s)),'°'];
    legEl2 = ['org \phi_',num2str(s),'=',num2str(elwDeg(s)),'°'];
    l8 = legend(legAz1, legEl1, legAz2, legEl2);
    set(l8, 'Location', 'EastOutside')
    title1 = ['Angle est. of ',num2str(s),'. source'];
    title(title1,'fontsize',titFontSize );
    xlabel('Angle (°)','fontsize',labFontSize );
    ylabel('Probability (-)','fontsize',labFontSize );
    set(hFig8, 'Position', [ 70 10 700 180]);

    hFig9(s) = figure;
    hist(azNsv(s, efBW(1):efBW(2)),360)
    title('distribution of ...')
    xlabel('magnitude [-]'); ylabel('[frequency]')
    end
end
%% PLOT separated time data only
if plotAnMSTMSSo(3)==1  % plot sep. in time?
% PLOT sep. time data : all microphones
% for s = 1:sigNum
%     tSep = [(1:1:size(p1org,2))/Fs];
%     hFig20 = figure(101);
%     set(hFig20, 'Position', [ 70 10 800 675]);
%     
%     hF1(s) = subplot(micNum+1,sigNum,1+micNum*(s));
%     plot(tSep, p1(s,:));
%     title(['u_1, s_', num2str(s)])
%     ylabel('Magnitude (-)','fontsize',10)    
%     xlabel('Time (s)','fontsize',10)
%     
%     hF2(s) = subplot(micNum+1,sigNum,2+micNum*(s));
%     plot(tSep, p2(s,:));
%     title(['u_2, s_', num2str(s)])
%     xlabel('Time (s)','fontsize',10)
%     
% 
%     hF3(s) = subplot(micNum+1,sigNum,3+micNum*(s));
%     plot(tSep, p3(s,:));
%     title(['u_3, s_', num2str(s)])
%     xlabel('Time (s)','fontsize',10)
%     
%     
%     hF4(s) = subplot(micNum+1,sigNum,4+micNum*(s));
%     plot(tSep, p4(s,:));
%     title(['u_4, s_', num2str(s)])
%     xlabel('Time (s)','fontsize',10)
%     
%     %     set(hFig301, 'Position', [ 70 10 800 120]);
%     hF5(s) = subplot(micNum+1,sigNum,s);
%     plot(tSep, x1234(s,1:size(tSep,2)));
%     title(['convolved signal: x_', num2str(s)])
%     xlabel('Time (s)','fontsize',10)
%     ylabel('Magnitude (-)','fontsize',10)    
% 
%     hF = [hF1, hF2, hF3, hF4, hF5];
%     axis(hF, [0 1 -2 2])
% end
% PLOT sep. time data : one microphone
    
    tSep = [(1:1:size(p1org,2))/Fs];
    hFig21 = figure(101);
    set(hFig21, 'Position', [ 70 10 500 675]);
    colXt = 'k';
    scConT = 1;
    pOff = 0;
%     hF5 = subplot(5,1,1);
hF5 = subplot(3,1,1);
    plot(tSep, nor(x1234(s,1:size(tSep,2))), colXt, 'LineWidth', tW);
    title(['Convolved sound signal p_1'],'fontsize',titFontSize)
%     xlabel('Time (s)','fontsize',labFontSize)
%     ylabel('Normalized magnitude (-)','fontsize',labFontSize)      
        
%     hF1 = subplot(5,2,4);
    hF1 = subplot(3,4,6);

    plot(tSep, nor(p1(1,:)), colXt,'LineWidth', tW);
    title(['Separated sound source v_1'],'fontsize',titFontSize)
    ylabel('Normalized magnitude (-)','fontsize',labFontSize)  
%     xlabel('Time (s)','fontsize',labFontSize)
    
%     hF6 = subplot(5,2,3);
    hF6 = subplot(3,4,5);
    plot(tSep, nor(p1sg(1,pOff+1:pOff+Fs))*scConT,colXt, 'LineWidth', tW);
    title(['Original sound source s_1'],'fontsize',titFontSize)
%     ylabel('Magnitude (-)','fontsize',labFontSize )    
%     xlabel('Time (s)','fontsize',labFontSize )
    
%     hF2 = subplot(5,2,6);
        hF2 = subplot(3,4,8);
    plot(tSep, nor(p1(2,:)),colXt, 'LineWidth', tW);
    title(['Separated sound  source v_2'],'fontsize',titFontSize)
%     xlabel('Time (s)','fontsize',labFontSize )
%     ylabel('Magnitude (-)','fontsize',labFontSize )    

%     hF7 = subplot(5,2,5);
        hF7 = subplot(3,4,7);
    plot(tSep, nor(p1sg(2,pOff+1:pOff+Fs))*scConT,colXt, 'LineWidth', tW);
    title(['Original sound source s_2'],'fontsize',titFontSize)
    ylabel('Magnitude (-)','fontsize',labFontSize )    
%     xlabel('Time (s)','fontsize',labFontSize )
%     ylabel('Normalized magnitude (-)','fontsize',labFontSize)  

%     hF3 = subplot(5,2,8);
    hF3 = subplot(3,4,10);
    plot(tSep, nor(p1(3,:)),colXt, 'LineWidth', tW);
    title(['Separated sound source v_3'],'fontsize',titFontSize)
    xlabel('Time (s)','fontsize',labFontSize )
%     ylabel('Magnitude (-)','fontsize',labFontSize )    

%     hF8 = subplot(5,2,7);
    hF8 = subplot(3,4,9);
    plot(tSep, nor(p1sg(3,pOff+1:pOff+Fs))*scConT,colXt, 'LineWidth', tW);
    title(['Origianl sound source s_3'],'fontsize',titFontSize)
%     ylabel('Magnitude (-)','fontsize',labFontSize )    
    xlabel('Time (s)','fontsize',labFontSize )
    
%     hF4 = subplot(5,2,10);
    hF4 = subplot(3,4,12);
    plot(tSep, nor(p1(4,:)),colXt, 'LineWidth', tW);
    title(['Separated sound source v_4'],'fontsize',titFontSize)
    xlabel('Time (s)','fontsize',labFontSize)
%     ylabel('Magnitude (-)','fontsize',labFontSize )
    
%     hF9 = subplot(5,2,9);
    hF9 = subplot(3,4,11);
    plot(tSep, nor(p1sg(4,pOff+1:pOff+Fs))*scConT,colXt, 'LineWidth', tW);
    title(['Original sound source s_4'],'fontsize',titFontSize)
%     ylabel('Magnitude (-)','fontsize',labFontSize )    
    xlabel('Time (s)','fontsize',labFontSize )
%     ylabel('Magnitude (-)','fontsize',labFontSize )    

    
    hF = [hF1, hF2, hF3, hF4, hF5, hF6, hF7, hF8, hF9];
    set([hF1, hF2, hF5, hF6, hF7],'Xticklabel',[])   
    if simul ==1
        axis(hF, [0 1 -10.5 10.5])
    elseif simul == 0;
        axis(hF, [0 1 -1.1 1.1])
    end
end
%% save plot play ? (time/freq. data)
% save all separated sources
if spscOrg(1)==1;
    save('matSep\SEP030614.mat', 'p1', 'p2', 'p3', 'p4', 'p4b', 'p5', 'p6', 'p7')
end
% plot original sources
if spscOrg(2) == 1;
    t = t(1:Fs);
    % plot signals on one plot
    hFig1 = figure;
    hFig2 = figure;
    for s = 1:sigNum
        figure(hFig1);
        plot(t, p1sg(s,:), 'LineWidth',tfW, 'color', colorLine(s,:)); hold on
        figure(hFig2);
        plot(f, Xspl(s,:), 'LineWidth',tfW, 'color', colorLine(s,:)); hold on
        leg{s} = ['s', num2str(s)];
        if s == sigNum
            figure(hFig1);
            plot(t, bk1s, 'LineWidth',tfW, 'color', colorLine(20,:));
            htTt = title(['Time course of ', num2str(sigNum),' unmixed signals'],'fontsize',titFontSize);
            xlabel('Time (s)','fontsize',labFontSize );
            ylabel('Magnitude (-)','fontsize',labFontSize );
            hLeg1 = legend(leg, 'bkn');
            set(hLeg1, 'location', 'northeast')
            [ytMin, ytMax, xtMin, xtMax] = axUni(t,p1sg, 1.1, 1);
            axis([ytMin, ytMax, xtMin, xtMax])
            set(hFig1, 'Position', [ 70 380 500 300]);
            
            figure(hFig2);
            hFbk = plot(f, BKspl, 'LineWidth',tfW, 'color', colorLine(20,:));
            htTf = title(['Frequency spectrum of separated ', num2str(sigNum),'microphones'],'fontsize',titFontSize);
            xlabel('Frequency (Hz)','fontsize',labFontSize);
            ylabel('Magnitude (dB SPL)','fontsize',labFontSize);            
            hLeg2 = legend(leg, 'bkn');
            uistack(hFbk,'bottom');
%             set(hLeg2, 'location', 'southwest')
            [yfftMin, yfftMax, xfftMin, xfftMax] = axUni(f,Xspl, 1, 0);
            axis([yAx xfftMin, xfftMax])
            set(hFig2, 'Position', [ 570 380 500 300]);
            x = [0 0 efBW(1) efBW(1) 0]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
            h1 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
            x = [efBW(2) efBW(2) efBW(2)+efBW(1) efBW(2)+efBW(1) efBW(2)]; y1 = [-efBW(1) efBW(1) efBW(1) -efBW(1) -efBW(1)]; % define edges of area
            h2 = patch(x, y1, -1 * ones(size(x)), [0.3 0.54 0.58], 'LineStyle', 'none');
%     set(htTt,'String','')
%     set(htTf,'String','')
%     picName1 = ['timeUnmix'];
%     picName2 = ['fftUnmix'];
%     picForm = 'eps';    % jpeg
%     picPath = 'refresh270514\figs\';
%     hgexport(hFig1, [picPath, picName1], hgexport('factorystyle'), 'Format', picForm);
%     hgexport(hFig2, [picPath, picName2], hgexport('factorystyle'), 'Format', picForm);        
        end
    end
    % plot signals separately
    for s = 1:sigNum %
        leg{s} = ['p', num2str(s)];
        hFig3(s) = figure;
        plot(t, p1sg(s,:), 'color', colorLine(s,:));
        title(['Time course of mic ', num2str(s)],'fontsize',titFontSize);
        xlabel('Time (s)','fontsize',labFontSize);
        ylabel('Magnitude (-)','fontsize',labFontSize);
        axis([ytMin, ytMax, xtMin, xtMax])
        hLsep1 = legend(leg{s});
%         set(hLeg2, 'location', 'southwest')
        set(hFig3(s), 'Position', [ 70 10 400 220]);
        
        hFig4(s) = figure;
        [X(s,:), f, t, hAn1] = anAudio(p1sg(s,:), Fs, 0);
        Xspl(s,:) = abs(20*log10((X(s,:))/2e-5));
        plot(f, Xspl(s,:), 'LineWidth',1, 'color', colorLine(s,:));
        title(['Frequency spectrum of mic ', num2str(s)],'fontsize',titFontSize);
        xlabel('Frequency (Hz)','fontsize',labFontSize);
        ylabel('Magnitude (dB)','fontsize',labFontSize);
        axis([yAx, xfftMin, xfftMax])
        hLsep2 = legend(leg{s});
%         set(hLeg2, 'location', 'southwest')
        set(hFig4(s), 'Position', [ 480 10 400 220]);        

        % Plot STFT
        B = 512;        % Block size
        Ra = 16;        % Number of samples between blocks
        Xst(s,:,:) = stftM(p1sg(s,:), B, Ra);  % Compute STFT
        fst = linspace(0, Fs/2, size(Xst(s,:,:),3));
        tst = linspace(0, length(p1sg(s,:))/Fs, size(Xst(s,:,:),2));
%         hFst(s) = figure; % in case spectrum alone
%         imagesc(tst, fst, log(abs(squeeze(Xst(s,:,:)))),[-8 6]); % 
%         axis xy
%         % ylim([0 20e3])
%         xlabel('Time (s)');
%         ylabel('Frequency (Hz)');
%         colorbar
        hFig5(s) = figure;
        hS1 = subplot(2,1,1);
        imagesc(tst, fst, log(abs(squeeze(Xst(s,:,:)))),[-8 6]);   % [-9.9 6]
        axis xy
        ylabel('Frequency (Hz)');
        set(hS1,'Xticklabel',[])        
        colorbar
        hS2 = subplot(2,1,2);
        plot(t(1:length(p1sg(s,:))), p1sg(s,:), 'b');
        xlabel('Time (s)');
        ylabel('Magnitude (-)');
        axis([ytMin, ytMax, xtMin, xtMax])
        % axis([0 Tmax, -max(abs(x1))-(max(abs(x1))*0.05) max(abs(x1))+(max(abs(x1))*0.05)])
        set(hFig5(s), 'Position', [ 920 10 400 250]);
        % get( axes2, 'XTick',  [0 0] );    % get rid of axis ticks
        ax1=get(hS1,'Position'); % [0.1300, 0.5838, 0.6512, 0.3412]
        ax2=get(hS2,'Position'); % [0.1300, 0.1100, 0.7750, 0.3412]
        % [shifts sideways, shifts up-down, stretches rightwards, stretches up]
        set(hS1,'Position',[ax1(1), ax1(2)-0.24, ax1(3)+0.124, ax1(4)+0.30]);
        set(hS2,'Position',[ax2(1), ax2(2), ax2(3), ax2(4)-0.1]);
%     picName = ['spectSig', num2str(sigN)];
%     picForm = 'eps';    % jpeg
%     picPath = 'Sound_source_separation\BSS\convolutive\sepSave\fig\';
%     hgexport(hFig, [picPath, picName], hgexport('factorystyle'), 'Format', picForm);
    end
end
% play all 
if spscOrg(3) == 1; 
    psgMix = zeros(size(p1sg(s,:)));
    for s = 1:sigNum
        sound(p1sg(s,:), Fs)
        pause(size(p1sg(s,:), 2)/Fs)
    end
    sound(p1org(1,:), Fs)
end
end
%% SAVE figures
% % save data?
% if spscOrg(1) == 1
%     save(sepData, 'x1234', 'p1sg', 'Fs');
% end
if saveEps == 1;
    % save angle of mix, single, time, multi, single multi detail, single only
    if plotAnMSTMSSo(1) == 1;   % mix angle
        set(hT30,'String','')
        set(hT40,'String','')        
        picName = 'mixAnglesLab';
        picForm = 'eps';    % jpeg
%         picPath1 = 'd:\Dropbox\_Uni-KAIST\3D_Intensity_array\docs\2014\3-5.LATEX\figures\';
        picPath = 'figs\';
        hgexport(hFig30, [picPath, picName], hgexport('factorystyle'), 'Format', picForm);
    end
    if saveFig == 1
    for s = 1:sigNum
%         set(hT70(s),'String','')
%         set(hT80(s),'String','')
        picName = ['singleAngleSLab',num2str(s)];
        picForm = 'fig';    % jpeg
        picPath = 'figs\';
        saveas(hFig70(s), [picPath, picName], picForm);
    end
    end
    if plotAnMSTMSSo(2) == 1;   %  single angle
        for s = 1:sigNum
        set(hT70(s),'String','')
        set(hT80(s),'String','')
        picName = ['singleAngleSLab',num2str(s)];
        picForm = 'eps';    % jpeg
        picPath = 'figs\';
        hgexport(hFig70(s), [picPath, picName], hgexport('factorystyle'), 'Format', picForm);
        end
    end
    if plotAnMSTMSSo(3) == 1;   % sep time data
%         set(hT20,'String','')
        picName = ['timeSepLab'];
        picForm = 'eps';    % jpeg
        picPath = 'figs\';
        hgexport(hFig20, [picPath, picName], hgexport('factorystyle'), 'Format', picForm);
    end
    if plotAnMSTMSSo(4) == 1;   % multiple separation angle
        set(hT10,'String','')
        set(hT20,'String','')
        picName = ['multiAnglesLab'];
        picForm = 'eps';    % jpeg
        picPath = 'figs\';
        hgexport(hFig10, [picPath, picName], hgexport('factorystyle'), 'Format', picForm);
    end
    if plotAnMSTMSSo(5) == 1;   % single multi angle detail
        for s = 1:sigNum
        set(hT50(s),'String','')
        set(hT60(s),'String','')
        picName = ['singleAngleSLab',num2str(s),'detail'];
        picForm = 'eps';    % jpeg
        picPath = 'figs\';
        hgexport(hFig50(s), [picPath, picName], hgexport('factorystyle'), 'Format', picForm);
        end
    end
end
%% table
if saveTab == 1
sourceName = {'$s_1$';'$s_2$';'$s_3$';'$s_4$'};
% make column vectors for variables to be put to tab
SNR_s = roundn(SNRband,-1)'; azDeg = roundn(azwDeg,-1)'; az_arit = roundn(azNm_arit,-1)'; az_err = roundn(azNm_err,-1)'; az_std = roundn(azNm_std,-1)'; 
elDeg = roundn(elwDeg,-1)'; el_arit = roundn(elNm_arit,-1)'; el_err = roundn(elNm_err,-1)'; el_std = roundn(elNm_std,-1)';
% create table
tabCvs = table(sourceName, SNR_s, azDeg, az_arit, az_err, az_std, elDeg, el_arit, el_err, el_std)%,  'VariableNames', varName') % 'RowNames', sourceName,
% print table to excel file
xlsName = ['tabs/',tabName, '.xlsx'];
if simul == 0
    [stat,popW] = xlswrite(xlsName,{'experimental data'},1,'B1');
    writetable(tabCvs,xlsName,'Sheet',1,'Range','B3:K7')
    if stat == 0
        display(['error: ', popW]);
    end
elseif simul ==1
    [stat,popW] = xlswrite(xlsName,'simulation data',2,'B1');
    writetable(tabCvs,xlsName,'Sheet',2,'Range','B3:K7')
end
% warning('off','MATLAB:xlswrite:AddSheet') % if new spredsheet created
end

% Functions
function yNor = nor(x)
    % normalize input by its absolute maximum
    yNor = x./max(abs(x));
end

% EOF