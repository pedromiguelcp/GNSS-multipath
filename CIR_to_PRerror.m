clear all;
close all;
clc;

% Examples:
chip_to_seconds=0.001/1023;
Channel1.Delay=0.*chip_to_seconds;
Channel1.Power=0;
Channel1.Phase=0;
RangeError=calcPseudoRangeError(Channel1); % multipath free
disp(RangeError);

Channel2.Delay=[0 0.1 0.23].*chip_to_seconds;
Channel2.Power=[1.01 -3 -4.2];
Channel2.Phase=[0.10 2.66 1.4];
RangeError=calcPseudoRangeError(Channel2); % LOS + 2 echos
disp(RangeError);

Channel3.Delay=[0 0.1 0.12 0.2 0.203].*chip_to_seconds;
Channel3.Power=[-1.07 -3 -3.2 -4 -4.01];
Channel3.Phase=[-2.2 0.25 2 0.5 0.09];
RangeError=calcPseudoRangeError(Channel3); % LOS + 4 echos
disp(RangeError);




function [PRerror] = calcPseudoRangeError(CIR)
    % Function to compute the pseudo range error from the Channel Impulse Response
    %
    % [PRerror] = calcPseudoRangeError(CIR)
    %
    %   Inputs:
    %       CIR - Channel Impulse Response
    %           CIR.Delay - CIR delays (seconds - relative to LOS)
    %           CIR.Power - CIR power  (dB)
    %           CIR.Phase - CIR phase  (radians)
    %
    %   Outputs:
    %       PRerror - Pseudo Range Error (meters)


    
    % Correlation vars
    c = 299792458; % m/s
    fs = 26000000; % Sample rate (Hz)
    codeFreq = 1.023e6; % C/A code rate (chips/s)
    codeFreqBasis = codeFreq;
    codelength = 1023;
    numSample = round(codelength/(codeFreq/fs));
    Correlator_spacing = 0.5; % wide correlator
    Spacing = [-Correlator_spacing 0 Correlator_spacing];
    tau2code = 0.3749245;
    tau1code = 0.007030542258775;
    remChip = 0;
    code_outputLast = 0;
    DLLdiscriLast = 0;
    codeError = 0.001;
    INCode   = 0;% LOS signal

    % CIR unit convertion
    ms = 0.001;
    seconds_to_chips=codelength/ms;
    CH.Delay=CIR.Delay.*seconds_to_chips;
    CH.Power=db2pow(CIR.Power);
    CH.Phase=CIR.Phase;

    % Generate the C/A code
    Code = generateCAcode(1); % satellite 1
    Code = [Code(end-1) Code(end) Code Code(1) Code(2) Code(3)]; % replica slidding
    
    % Interference injection
    for idx = 1: size(CH.Delay, 2)
        multipath = Spacing(2) + CH.Delay(idx) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(2) + CH.Delay(idx));
        NLOS_to_LOS_phase_ratio=cos(CH.Phase(1)-CH.Phase(idx));
        INCode = INCode + Code(ceil(multipath) + 2).*(CH.Power(idx)*NLOS_to_LOS_phase_ratio);
    end
    
    while (abs(codeError) >= 0.001) % convergence threshold
        % Correlator engine
        codePhaseStep   = codeFreq/fs;
        t_CodeEarly     = (Spacing(1) + remChip) : codePhaseStep : ((numSample -1) * codePhaseStep + Spacing(1) + remChip);
        t_CodePrompt    = (Spacing(2) + remChip);
        t_CodeLate      = (Spacing(3) + remChip) : codePhaseStep : ((numSample -1) * codePhaseStep + Spacing(3) + remChip);
        CodeEarly       = Code(ceil(t_CodeEarly) + 2);
        CodeLate        = Code(ceil(t_CodeLate) + 2);
        E_i             = sum(CodeEarly  .*INCode);
        L_i             = sum(CodeLate   .*INCode);
        remChip         = (((numSample -1) * codePhaseStep + t_CodePrompt) + codePhaseStep) - codeFreqBasis*ms;

        % Discriminator
        codeError       = 0.5 * (abs(E_i)-abs(L_i))/(abs(E_i)+abs(L_i));
        codeNco         = code_outputLast + (tau2code/tau1code)*(codeError - DLLdiscriLast) + codeError* (ms/tau1code);
        DLLdiscriLast   = codeError;
        code_outputLast = codeNco;
        codeFreq        = codeFreqBasis - codeNco;
    end

    % Pseudo range error based on correlators final state of convergence - more accurate
    PRerror=((Spacing(2) + remChip)-CH.Delay(1))*(ms/codelength)*c;
end