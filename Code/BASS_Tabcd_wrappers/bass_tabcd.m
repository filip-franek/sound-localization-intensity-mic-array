%% bass_tabcd.m Blind Audio Sound Source Separation (BASS) function
%   
%	This wrapper returns separated source from audio mixture alowing user to pick a BASS method originally available
%	at tddeconv [1]. Thus making it possible to use the BASS algorithm in an mfile.
%     
%	Input variables
%   	x               :   input data of size(x)=[m N], m microphones, N samples
%       Fs              :   sampling frequency
%       ICAmet          :	method option:
%           bgl         :   Blog Gausian Likelihood
%           efica       :   Eficient Fast ICA
%           extefica    :   Block EFICA 
%           bwasobi     :   actually sort of BARBI Block-AutoRegressive Blind Identification
%           numBands    :   number of bands
%                       :   1 band......1
%                       :   2 bands.....2   
%                       :   4 bands.....3
%                       :   8 bands.....4
%                       :   16 bands....5
%       similStr        :   similarity comparison
%           'pro'       :   projections
%           'gcc'       :   GCC-PHAT
%           'coh'       :   Coherence
%       weiTypStr       :   weighting type
%           'norm'      :   'Normal'
%           'toep'      :   'Toeplitz'
%           'bin'       :   'Binary'
%           'fuzz'      :   'Fuzzy'
%           'fCon'      :   'F. Constrained'
%           'fBin'      :   'F.Binary'
%       clust           :   clusterring
%           'hclus'     :   hierarchical clustering
%           'rfcm'      :   clustering (fuzzy?)
%       weiTypStr       :   weighting type string
%       weiVal          :   weighting value
%       xName           :   name of processing input variable
% 
%   Output variables
%       y               :   original multichannel recording data
%       est             :   estimated seperated sources - est(i,:,j) ith microphone, jth source
%       data            :   parameters used
%       Xdis            :   ?
% 
% 
%   Usage: 
%          [y, est, data, Xdis] = bass_tabcd(x, 44100, 'bgl', 20,   'pro',       16 ,      1, 'hclus', 'norm',   1   , 'x');
% function [y, est, data, Xdis] = bass_tabcd(x, Fs, ICAmet, fiLen, similStr, subBandsR, muPar, clust, weiTypStr, weiVal, xName);
%                                            1  2     3       4      5      6         7      8      9       10      11  

% 
% Date: 19/01/2015
% Author: Filip Franek, AudibleBits, filip@audiblebits.com
% 
% Code modified from function tabcd obtained from https://asap.ite.tul.cz/downloads/t-abcd-a-time-domain-method-for-
%                       blind-audio-source-separation-based-on-a-complete-ica-decomposition-of-an-observation-space/
%
%  References:
%   [1] Koldovsky, Zbynek, and Petr Tichavsky. "Time-domain blind separation of audio sources on the basis of a complete 
%       ICA decomposition of an observation space." IEEE transactions on audio, speech, and language processing 19, 
%       no. 2 (2010): 406-416.
% 
% Original work Copyright (C) 2016 Zbynek Koldovsky and Petr Tichavsky
% Modified Work Copyright (C) 2020 Filip Franek
%%
function [y, est, data, Xdis] = bass_tabcd(x, Fs, ICAmet, fiLen, similStr, subBandsR, muPar, clust, weiTypStr, weiVal, xName);

    clear y est

    if size(x,1)>size(x,2) % check if row vector
        x=x';
    end

    if nargin==0
       error('Missing input signals!');
    end
    switch ICAmet
        case 'efica'
            ICAnum=1;
        case 'bgl'
            ICAnum=2;
        case'extefica'
            ICAnum=3;
        case 'bwasobi'
            ICAnum=4;
        otherwise fprintf('!!!!! undefined ICA method !!!!!!\n')
    end
    switch similStr
        case 'pro'
            simil = 1;
        case 'coh'
            simil = 2;
        case 'gcc'
            simil = 3;
        otherwise fprintf('wrong similarity method\n')    
    end
    switch subBandsR
        case 1
            subBands = 1;
        case 2
            subBands = 2;
        case 4
            subBands = 3;
        case 8
            subBands = 4;
        case 16
            subBands = 5;
        otherwise fprintf('undefined number of bands\n')
    end
    switch weiTypStr
        case 'norm'
            weiTyp= 1;
        case 'toep'
            weiTyp= 2;
        case 'bin'
            weiTyp= 3;
        case 'fuzz'
            weiTyp= 4;
        case 'fCon'
            weiTyp= 5;
        case 'fBin'
            weiTyp= 6;
        otherwise fprintf('undefined number of bands\n')        
    end

    data.x=x;                   % 1
    data.fs = Fs;               % 2
    data.method = ICAmet;       % 3
    data.filterlength = fiLen;  % 4
    data.nsubbands=subBands;    % 5
    data.mu = muPar;            % 6
    data.clustering=clust;      % 7
    data.similarity = simil;    % 8
    data.wtype = weiTyp;        % 9
    data.weighting = weiVal;    % 10
    data.lengthofplot=floor(length(x));
    data.beginofplot=1;
    data.lengthofICA=length(x);
    data.offsetICA=1;
    data.iweights=0;
    data.filterlengtheffective=(1+0.4*abs(data.mu-1)*log10(data.filterlength))*data.filterlength/data.mu;
    
      str = sprintf('T-ABCD BBS, method: %s, signal:%s',ICAmet,xName);
      figNumber=figure( ...
         'Name',str, ...
         'RendererMode','manual',...
         'Position',[100,100,600,600],...
         'Visible','off');

          param.lengthofICA=data.lengthofICA;
          param.offsetICA=data.offsetICA;
          param.method=data.method;
          param.weighting=data.weighting;
          param.wtype=data.wtype;
          param.clustering=data.clustering;
          param.mu=data.mu;
          param.useiweights=data.iweights;
          param.similarity=data.similarity;


          if data.nsubbands==1
              [data.shat data.microphones We Xdis]=bass_deconv(data.x,data.filterlength-1,param);
          else % subband separation
              load gfilter          
              nbands=2^(data.nsubbands-1);
              N1=ceil(size(data.x,2)/nbands);
              U=zeros(nbands,N1,size(data.x,1));
              for i=1:size(data.x,1)
                  U(:,:,i)=treeanalysis(g,data.x(i,:),data.nsubbands-1);
              end
              SBres=zeros(size(data.x,1),size(U,2),size(data.x,1),size(U,1));
              for i=1:size(U,1)
                  fprintf('Band number %d ------------\n',i);
                  % TODO: Add condition so that number of sources is equal to number of microphones
                  param.offsetICA=ceil(data.offsetICA/nbands);
                  param.lengthofICA=ceil(data.lengthofICA/nbands);
                  [~, SBres(:,:,:,i), We, Xdis]=bass_deconv(squeeze(U(i,:,:))',data.filterlength-1,param);              
              end
              SBres=SBpermute3(SBres);
              if isfield(data,'microphones')
                  data=rmfield(data, 'microphones');
              end
              for m=1:size(data.x,1)
                for n=1:size(data.x,1)
                 data.microphones(m,:,n)=treesynthesis(g,squeeze(SBres(m,:,n,:))');
                end
              end
              data.shat=squeeze(data.microphones(1,:,:))';
          end

          est = data.microphones;
          y = data.shat;
    end