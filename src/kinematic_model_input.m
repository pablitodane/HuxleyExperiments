function [stim,lmtc,lmtcd]=kinematic_model_input(t,parms)
%% Original code - Stimulation
% if t<=.1
%     stim=parms.gamma0; % we start in steady state
% elseif t>.1 && t<2 % full activation
%     stim=1; 
% else
%     stim=0.2; % relaxation to low value
% end
% lmtc0=parms.lmtc0;
% lceopt=parms.lceopt;

%% Modified code - Stimulation
if t<=0.01
    stim=parms.gamma0; % we start in steady state
elseif t>0.01 && t<0.85 % full activation
    stim=1; 
else
    stim=0.0; % relaxation to zero
end
lmtc0=parms.lmtc0;
lceopt=parms.lceopt;

% mtc length
% isometric at lmtc0

%% Original code - Kinematic input
% if t<=0.3
%     lmtc=lmtc0;
%     lmtcd=0;
% else
% 
% lmtc=lmtc0+0.2*lceopt*sin(2*pi*t);
% lmtcd=0.2*2*pi*lceopt*cos(2*pi*t);
% end

%% Modified code - Kinematic input
if t<=0.45
    lmtc=lmtc0;
    lmtcd=0;
elseif t>0.45 && t<=0.95
    lmtc=lmtc0 - 0.001*0.45 + 0.001*t;
    lmtcd=0.001;
else
    lmtc=lmtc0 - 0.001*0.45 + 0.001*1;
    lmtcd=0;
end

end