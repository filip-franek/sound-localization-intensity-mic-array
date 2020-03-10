%% multi_si_bass.m Multiple Sound Source Localization Using 3D Sound Intensity Method and BASS
%
%	Sound sources location estimation based on 3D sound intensity measurement at presence of 4 sound sources 
%   simultaneously playing. Blind Audio Source Separation (BASS) is used to separate the sound sources.
%   Sound recordings are taken in an anechoic chamber.
%   
%   clone github repository and run this file to process measured data files and return plots
%
%   variables
%   	overWrite   :   overwrite BASS results or load from mat file [true/false]
% 
% Date: 10/03/2020
% Author: Filip Franek, AudibleBits, filip@audiblebits.com
%
%  References:
%   [1] Fahy, Frank J., and Vincent Salmon. "Sound intensity." (1990): 2044-2045.
%   [2] Koldovsky, 
%   [3] Franek, Filip. "Localization of Incoherent Sound Sources by
%       Three-dimensional Intensity Array." (2014)
%
% Original work Copyright (C) 2020 Zbynek Koldovsky and Petr Tichavsky
% Modified Work Copyright (C) 2020 Filip Franek
%%

clear all; clc; close all;

% Add path to subfunctions
addpath('BASS_Tabcd_wrappers')
addpath('BASS_Tabcd')

% Add script's folder to Matlab path
root = mfilename('fullpath');
[rootDir,name,ext] = fileparts(root);
cd(rootDir);

% Measurement data file
mixFile = 'sp1234mix_0dB';  % 4 microphone in tetrahedral config recorded 4 speakers simultaneously
roomExp = 'anCham';         % Experiment environment

azwDeg = [-12.2, 54.4, 29.2, -34.8 ]; % True azimuth angle of sound sources [1 2 3 4]
elwDeg = [11.8, -15.3, -10.7, -45.5]; % True elevation angle

% Measured background noise in measurement location
bkFile = 'bkNoise';

d=50*0.001;         % distance between microphones in tetrahedral configuration
Fs = 2^13;          % audio sampling freq. (2^13=8192 Hz)

% Make speech always first for further comparison
sigTyp = {'Speech1', 'Normal', 'Speech3', 'Speech4'};   % 'Normal', 'Laplac', 'Uniform'
% Take section of audio to process
secL = 1.5;     % secL always larger then secLx
secLx = 1;
% if secL=<secLx
samLen = Fs*secL;
samLenx = Fs*secLx;
sigNum = length(sigTyp);    % number of sources
% Constants
c = 343;            % Speed of sound
rho=1.21;           % Density of air
df = 2^(0);         % frequency resolution in Hz (2^n, n = (-3:1:+6))
bw = 100;            % bandwidth over which averaged is taken
nfft = Fs/df;       % number of FFT points
%% Plot options
% M - multi, S - sigle, X - mix, d - detail, 0 - original angle, 
angMdW = 1.5;    angMdWo = 2;    angMdS = '-';     % W - lineWidth, S - lineStyle
angMW  = 1.0;    angMWo  = 0.5;  angMS  = '-'; angMoS = '--';      
angXW  = 1.5;    angXWo  = 2;    angXS  = '.';
sgRan = 90;             % angle range for single unmixed source
sgRanS = 90;            % angle range for single separated sources : ylim +-X degrees
mixRan = 90;            % angle range of mixed not separated signals
xOver = 0;              % plotting x below and over effective bandwidth
efBW = [500, 1500];     % effective freq. bw. of intensity array
yAx = [efBW(1)-xOver, efBW(2)+xOver]; % xlimit for frequency
labFontSize = 10;
legFontSize = 8;
titFontSize = 10;
titFontWeight = 'bold';
legFontSizeX = 10;  % because mixed seems different to others
tfW = 1;	% time-frequency plot, width
tW = 1;

sepFold = 'matSep';

overWrite = 0;      % overwrite BASS results or load from mat file !!!!!!!!?

%% ++++++++++++MULTI main+++++++++++++
%%
% TODO: Seperate source signals fn
% TODO: Calculate SI fn
% TODO: Plots fns

% TODO: Load microphone data fn
% input:
%           matPath
%           onset
%           lenS
%% LOAD Multi Data
% if simul == 0;
matFold = '../Data/matOrg/';
matPath = [matFold,mixFile,'.mat'];
% structure of sp1234mix_0dB

load(matPath)

onset = 0;% in sec
lenS = 1;% in sec
phpTrun = php(onset*Fs+1:(onset+lenS)*Fs,:)';

% forming microphone group
x1234 = phpTrun(1:4,:);  x1234Name = varname(x1234);
x4567 = phpTrun(4:7,:);  x4567Name = varname(x4567);
x1234567 = phpTrun(1:7,:);  x1234567Name = varname(x1234567);
if overWrite == 1
    ICAmet = 'efica';   % bgl efica extefica bwasobi
    simil = 'gcc';      %1...projections 2...GCC-PHAT 3...Coherence
    subBandsR = 1;
    muPar = 1;
    clust = 'hclus';    % cluster: rfcm=fuzzy c-mean, hclus=hierarchical
    weiTyp = 'toep';	% varies-1:normal,2:toeplitz,3: binary...
    weiVal = 1;    
    fiLen = 20;         % filter length
    [y1, est1, data1, Xdis1] = bass_tabcd(x1234(:,onset+1:Fs+onset),Fs, ICAmet, fiLen, simil, subBandsR, muPar, clust, weiTyp, weiVal, x1234Name);
    [y2, est2, data2, Xdis2] = bass_tabcd(x4567(:,onset+1:Fs+onset),Fs, ICAmet, fiLen, simil, subBandsR, muPar, clust, weiTyp, weiVal, x4567Name);

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

    indC1 = [2, 4, 3, 1]; % [4, 1, 3, 2]
    indC2 = [3, 4, 1, 2]; % [3, 4, 1, 2]
    % FOR     'efica','gcc','hclus','norm'
    indC1 = [4, 2, 3, 1]; % [4, 2, 3, 1]
    indC2 = [1, 2, 4, 3]; % [1, 2, 4, 3]

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
    load('../Data/matSep\SEPanec270514.mat', 'p1', 'p2', 'p3', 'p4', 'p4b', 'p5', 'p6', 'p7')
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
    % module 1
    [G21 fout] = cpsd(p2(s,:),p1(s,:),nfft,0,nfft,Fs,'onesided');
    [G31 ~] = cpsd(p3(s,:),p1(s,:),nfft,0,nfft,Fs,'onesided');
    [G41 ~] = cpsd(p4(s,:),p1(s,:),nfft,0,nfft,Fs,'onesided');
    [G32 ~] = cpsd(p3(s,:),p2(s,:),nfft,0,nfft,Fs,'onesided');
    [G42 ~] = cpsd(p4(s,:),p2(s,:),nfft,0,nfft,Fs,'onesided');
    [G43 ~] = cpsd(p4(s,:),p3(s,:),nfft,0,nfft,Fs,'onesided');
    % module 2
    [G45 ~] = cpsd(p4b(s,:),p5(s,:),nfft,0,nfft,Fs,'onesided');
    [G46 ~] = cpsd(p4b(s,:),p6(s,:),nfft,0,nfft,Fs,'onesided');
    [G47 ~] = cpsd(p4b(s,:),p7(s,:),nfft,0,nfft,Fs,'onesided');
    [G65 ~] = cpsd(p6(s,:),p5(s,:),nfft,0,nfft,Fs,'onesided');
    [G75 ~] = cpsd(p7(s,:),p5(s,:),nfft,0,nfft,Fs,'onesided');
    [G76 ~] = cpsd(p7(s,:),p6(s,:),nfft,0,nfft,Fs,'onesided');    

    g21 = G21; g31 = G31; g32 = G32; g41 = G41; g42 = G42; g43 = G43;
    g45 = G45; g46 = G46; g47 = G47; g65 = G65; g75 = G75; g76 = G76;

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
    azNm(s,:)=atan2(IYm(s,:),IXm(s,:))*(180/pi);
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
    azNm1(s,:)=atan2(IYm1(s,:),IXm1(s,:))*(180/pi);
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
    azNm2(s,:)=atan2(IYm2(s,:),IXm2(s,:))*(180/pi);
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
    hT10 = title(title_azN1,'fontsize',titFontSize, 'fontweight', titFontWeight);
    xlabel('Frequency (Hz)','fontsize',labFontSize );
    ylabel('Azimuth angle (°)','fontsize',labFontSize );
    set(hFig10, 'Position', [ 70 10 800 620]);
    if s ==sigNum   % legend for 2 sources
        a = 'Ignoring extra legend entries.';
        MSGID ='MATLAB:legend:IgnoringExtraEntries';
        warning('off', MSGID)
        l1 = legend(legAz1{1},legAz2{1},legAz1{2},legAz2{2},legAz1{3},legAz2{3},legAz1{4},legAz2{4}, 'fontsize',legFontSize );
        
        set(l1,'Location','EastOutside');
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
    hT20 =title(title_azN1,'fontsize',titFontSize );
    xlabel('Frequency (Hz)','fontsize',labFontSize )
    ylabel('Elevation angle (°)','fontsize',labFontSize )
    if s ==sigNum   % legend for 2 sources
        l2 = legend(legEl1{1},legEl2{1},legEl1{2},legEl2{2},legEl1{3},legEl2{3},legEl1{4},legEl2{4}, 'fontsize',legFontSize);
        set(l2,'Location','EastOutside');
        
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

%% Plot histogram
hFig60 = figure(60);    % histogram multiple sources
histRng = 15;           % angle display range +-histRng
binRng = 0.1;           % angle resolution
hFig65 = figure(65);  % 3D histogram + 2D plot azimuth

for s = 1:sigNum;
    azHist(s,:) = (azwDeg(s)-180):binRng:(azwDeg(s)+180);       % set histogram vector
    elHist(s,:) = (elwDeg(s)-180):binRng:(elwDeg(s)+180);
    binCount(s,:) = histc(azNm(s,efBW(1):efBW(2)),azHist(s,:)); % compute histogram for azi, no plot
    binCount(s+4,:) = histc(elNm(s,efBW(1):efBW(2)),elHist(s,:));   % compute histogram for ele, no plot
end
[histVal histInd] = max(max(binCount));     % find maximum value of histogram
for s = 1:sigNum;
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
%     hSub66(s) = psubplot(2,4,s);
%     figure(hFig69)
%     hSub69(s) = subplot(1,4,s);
    binRng3d = 1;	% width for x and y bins
    data = [azNm(s,efBW(1):efBW(2));elNm(s,efBW(1):efBW(2))]';
    x = [azwDeg(s)-histRng:binRng3d:azwDeg(s)+histRng ];	% use same range for y bins, can be different
    y = [elwDeg(s)-histRng:binRng3d:elwDeg(s)+histRng ]; % range for x bins

    hSub67(s) = subplot(2,2,s);
    n = hist3(data, 'edges', {x,y}); % Extract histogram
    colormap(flipud(gray))
    hHex(s) = pcolor(x,y,n');
    xlabel(['Azimuth \theta_', num2str(s), '(°)'],'fontsize',labFontSize );
    ylabel(['Elevation \phi_', num2str(s), '(°)'],'fontsize',labFontSize );
    axis square
        set(gcf, 'renderer', 'opengl');
    set(get(gca,'child'), 'FaceColor', 'interp', 'EdgeColor','w', 'CDataMode', 'auto');

end
% % set(l61,'Location','EastOutside');
%     set(h,'FaceColor',[0 .5 .5],'EdgeColor','w')
set(hFig65, 'Position', [ 70 325 900 420]);
set(hFig60, 'Position', [ 70 10 900 420]);
sgtitle(hFig60, 'Histogram of separated sound sources', 'fontsize', titFontSize, 'fontweight', titFontWeight)
sgtitle(hFig65, '2D intensity map of separated sound sources', 'fontsize', titFontSize, 'fontweight', titFontWeight)

% Function to return variable name
function out = varname(var)
  out = inputname(1);
end