function [stim,lmtc,lmtcd]=kinematic_model_input(t,parms)
% This function can be modified to choose what kind of figure we want to have as an
% output. The options are:
% 1  - Original code from the tutorial
% 2 - Simulate the figure 4a
% 3 - Simulate the figure 4b
% 4 - Simulate the figure 4c
% 5 - Simulate the figure 4d

what_code = 3;

if what_code == 1
    
% Original code - Stimulation
    if t<=.1
        stim=parms.gamma0; % we start in steady state
    elseif t>.1 && t<2 % full activation
        stim=1; 
    else
        stim=0.2; % relaxation to low value
    end
lmtc0=parms.lmtc0;
lceopt=parms.lceopt;

% Original code - Kinematic input
    if t<=0.3
        lmtc=lmtc0;
        lmtcd=0;
    else
        lmtc=lmtc0+0.2*lceopt*sin(2*pi*t);
        lmtcd=0.2*2*pi*lceopt*cos(2*pi*t);
    end


elseif what_code == 2

% Figure 4a - Modified code - Stimulation
    if t<=0.01
        stim=parms.gamma0; % we start in steady state
    elseif t>0.01 && t<0.5 % full activation
        stim=1; 
    else
        stim=0.0; % relaxation to zero
    end
lmtc0=parms.lmtc0;
lceopt=parms.lceopt;    

% Figure 4a - Modified code - Kinematic input
    if t<=0.47
        lmtc=lmtc0;
        lmtcd=0;
    elseif t>0.47 && t<=0.48
        lmtc=lmtc0 + 0.01*0.47 - 0.01*t;
        lmtcd=-0.01;
    elseif t>0.48 && t<=0.83
        lmtc=lmtc0 + 0.01*0.47 - 0.01*0.48;
        lmtcd=0;
    elseif t>0.83 && t<=0.85
        lmtc=lmtc0 + 0.01*0.47 - 0.01*0.48 - 0.005*0.83 + 0.005*t;
        lmtcd=0.005;
    else
        lmtc=lmtc0 + 0.01*0.47 - 0.01*0.48 - 0.005*0.83 + 0.005*0.85;
        lmtcd=0;
    end


elseif what_code == 3

% Figure 4b - Modified code - Stimulation
    if t<=0.01
        stim=parms.gamma0; % we start in steady state
    elseif t>0.01 && t<0.657 % full activation
        stim=1;
    else
        stim=0.0; % relaxation to zero
    end
lmtc0=parms.lmtc0;
lceopt=parms.lceopt;

% Figure 4b - Modified code - Kinematic input
    if t<=0.497
        lmtc=lmtc0;
        lmtcd=0;
    elseif t>0.497 && t<=0.504
        lmtc=lmtc0 + 0.0143*0.497 - 0.0143*t;
        lmtcd=-0.0143;
    elseif t>0.504 && t<=1.002
        lmtc=lmtc0 + 0.0143*0.497 - 0.0143*0.504 + 0.00161*0.504 + -0.00161*t; 
        lmtcd=-0.00161;
    else
        lmtc=lmtc0 + 0.0143*0.497 - 0.0143*0.504 + 0.00161*0.504 + -0.00161*1.002;
        lmtcd=0;
    end


elseif what_code == 4

% Figure 4c - Modified code - Stimulation
    if t<=0.01
        stim=parms.gamma0; % we start in steady state
    elseif t>0.01 && t<0.85 % full activation
        stim=1; 
    else
        stim=0.0; % relaxation to zero
    end
lmtc0=parms.lmtc0;
lceopt=parms.lceopt;

% Figure 4c - Modified code - Kinematic input
    if t<=0.45
        lmtc=lmtc0;
        lmtcd=0;
    elseif t>0.45 && t<=0.95
        lmtc=lmtc0 - 0.001*0.45 + 0.001*t;
        lmtcd=0.001;
    else
        lmtc=lmtc0 - 0.001*0.45 + 0.001*0.95;
        lmtcd=0;
    end




elseif what_code == 5

% Figure 4d - Modified code - Stimulation
    if t<=0.01
        stim=parms.gamma0; % we start in steady state
    elseif t>0.01 && t<0.505 % full activation
        stim=1;
    else
        stim=0.0; % relaxation to zero
    end
lmtc0=parms.lmtc0;
lceopt=parms.lceopt;

% Figure 4d - Modified code - Kinematic input
    if t<=0.49026
        lmtc=lmtc0;
        lmtcd=0;
    elseif t>0.49026 && t<=0.8943
        lmtc=lmtc0 + 0.002475*0.49026 - 0.002475*t;
        lmtcd=-0.002475;
    else
        lmtc=lmtc0 + 0.002475*0.49026 - 0.002475*0.8943;
        lmtcd=0;
    end

end % end-if

end