function [file, signal, acq, track, solu, cmn] = initParameters()
%Purpose:
%   Parameter initialization
%Inputs: 
%	None
%Outputs:
%	file        - parameters related to the data file to be processed,a structure
%	signal      - parameters related to signals,a structure
%	acq         - parameters related to signal acquisition,a structure
%	track       - parameters related to signal tracking,a structure
%	solu        - parameters related to navigation solution,a structure
%	cmn         - parameters commmonly used,a structure
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.2
% 
% Written by B. XU and L. T. HSU

%% File parameters
file.fileName       = 'USRPN310_Scenario1_25M_I8Q8_100s'; 
file.fileRoute      = ['C:\Users\Pedro\Documents\phd\data\',file.fileName,'.dat']; 
file.skip        	= 0; % in unit of ms
solu.iniPos	= [41.451401319624345/180 * pi, -8.29113842707692/180 * pi, 50]; % Ground truth location
global ALPHA BETA 
ALPHA = [8.3819E-09  1.4901E-08 -5.9605E-08 -5.9605E-08];%       IONOSPHERIC CORR
BETA  = [8.3968E+04  1.6384E+04 -1.3107E+05 -6.5536E+04];%       IONOSPHERIC CORR
cmn.doy = 262; % Day of year
 

%% File parameters
file.fid           	= fopen(file.fileRoute,'r','ieee-le');
file.skiptimeVT     = 0; % skip time from the first measurement epoch of CT, in uint of msec
file.dataType       = 2;    %1:I; 2:IQ
file.dataPrecision  = 1;    %1:int8 or byte; 2; int16 

%% Signal parameters
signal.IF               = 0;%1580e6-1575.42e6; % unit: Hz 
signal.Fs               = 25e6;	% unit: Hz
signal.Fc               = 1575.42e6; % unit: Hz	
signal.codeFreqBasis	= 1.023e6; % unit: Hz 	
signal.ms               = 1e-3; % unit: s
signal.Sample           = ceil(signal.Fs*signal.ms);	
signal.codelength       = signal.codeFreqBasis * signal.ms;

%% Acquisition parameters
acq.prnList     = 1:32;	% PRN list
acq.freqStep    = 500;	% unit: Hz
acq.freqMin     = -10000;   % Minimum Doppler frequency
acq.freqNum     = 2*abs(acq.freqMin)/acq.freqStep+1;    % number of frequency bins
acq.L           = 10;   % number of ms to perform FFT

%% Tracking parameters
track.mode                  = 0;    % 0:conventional tracking; 1:vector tracking
track.CorrelatorSpacing  	= 0.05;  % unit: chip
track.DLLBW               	= 2;	% unit: Hz
track.DLLDamp           	= 0.707; 
track.DLLGain            	= 0.1;	
track.PLLBW              	= 15;
track.PLLDamp             	= 0.707;
track.PLLGain              	= 0.25; 	
track.msToProcessCT       	= 80000; % unit: ms 40000
track.msPosCT               = 10000; % unit: ms
track.msToProcessVT         = 2000; %track.msPosCT - file.skiptimeVT; %
track.pdi                   = 1; %


%% Navigation solution parameters
solu.navSolPeriod = 20; % unit: ms 
solu.mode  	= 0;    % 0:conventional LS/WLS; 1:conventionalKF; 2:VT


%% commonly used parameters
cmn.vtEnable  	= 0;%   % 0: disable vector tracking; 1:enable vector tracking
cmn.cSpeed      = 299792458;    % speed of light, [m/s]
global m2lat m2lon
m2lat = 1/110734;
m2lon = 1/103043; 
