clc, clear all

%%%% Let's state the main directories
rootDir = '/home/paolo/github/HuxleyExperiments/';
srcDir  = [rootDir,'src/'];
plotDir = [rootDir,'plots/'];
matDir  = [rootDir,'mat/'];
csvDir  = [rootDir,'csv/'];

%% Fig02
%%%% Fig02
fig2Dir = [matDir,'GraphComparison/fig2/'];

load([fig2Dir,'f_attach.mat']);
load([fig2Dir,'g_detach.mat']);
load([fig2Dir,'fx0.mat']);
load([fig2Dir,'gx0.mat']);
load([fig2Dir,'x0.mat']);

figure;plot(x0,fx0,x0,gx0,f_attach(:,1),f_attach(:,2),g_detach(:,1),g_detach(:,2))
title('Fig. 2 - rate parameter functions')
xlabel('x [h]')
ylabel('rate [Hz]')
xlim([-2 2])
legend('attachment - model','detachment - model','attachment - figure', 'detachment - figure')

%% Fig04
%%%% Fig04_A
fig4aDir = [matDir, 'GraphComparison/fig4/fig4a/'];
load([fig4aDir, 'fce.mat']);
load([fig4aDir, 't.mat']);
load([fig4aDir, 'fig4A.mat']);


figure;plot(t,[fce],fig4A(:,1),fig4A(:,2))
title('Fig 4a - forces')
legend('ce - model','ce - graph')
xlabel('Time [s]')
ylabel('Force [N]')

%%%% Fig04_B
fig4bDir = [matDir, 'GraphComparison/fig4/fig4b/'];
load([fig4bDir, 'fce.mat']);
load([fig4bDir, 't.mat']);
load([fig4bDir, 'fig4B.mat']);


figure;plot(t,[fce],fig4B(:,1),fig4B(:,2))
title('Fig 4b - forces')
legend('ce - model','ce - graph')
xlabel('Time [s]')
ylabel('Force [N]')

%%%% Fig04_C
fig4cDir = [matDir, 'GraphComparison/fig4/fig4c/'];
load([fig4cDir, 'fce.mat']);
load([fig4cDir, 't.mat']);
load([fig4cDir, 'fig4C.mat']);


figure;plot(t,[fce],fig4C(:,1),fig4C(:,2))
title('Fig 4c - forces')
legend('ce - model','ce - graph')
xlabel('Time [s]')
ylabel('Force [N]')

%%%% Fig04_D




%% Fig05
%%%% Fig05_A

%%%% Fig05_B

%%%% Fig05_C

%%%% Fig05_D


%% Fig06
%%%% Fig06_A

%%%% Fig06_B

%%%% Fig06_C



%% Fig07
%%%% Fig07_A

%%%% Fig07_B

%%%% Fig07_C

%%%% Fig07_D