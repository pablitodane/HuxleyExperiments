function [stim,lmtc,lmtcd]=kinematic_model_input(t,parms)
% This function can be modified to choose what kind of figure we want to have as an
% output. The options are:
% 1  - Original code from the tutorial
% 2 - Simulate the figure 4a
% 3 - Simulate the figure 4b
% 4 - Simulate the figure 4c
% 5 - Simulate the figure 4d

what_code = 5;

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

elseif what_code == 5;

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

% mtc length
% isometric at lmtc0

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