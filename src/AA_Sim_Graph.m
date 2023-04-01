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
%%%% Fig04_A_Data
fig4aDir = [matDir, 'GraphComparison/fig4/fig4a/'];
load([fig4aDir, 'fce.mat']);
fce4a = fce;
load([fig4aDir, 't.mat']);
t4a=t;
load([fig4aDir, 'fig4A.mat']);

%%%% Fig04_B_Data
fig4bDir = [matDir, 'GraphComparison/fig4/fig4b/'];
load([fig4bDir, 'fce.mat']);
fce4b = fce;
load([fig4bDir, 't.mat']);
t4b=t;
load([fig4bDir, 'fig4B.mat']);

%%%% Fig04_C_Data
fig4cDir = [matDir, 'GraphComparison/fig4/fig4c/'];
load([fig4cDir, 'fce.mat']);
fce4c=fce;
load([fig4cDir, 't.mat']);
t4c=t;
load([fig4cDir, 'fig4C.mat']);

%%%% Fig04_D_Data
fig4dDir = [matDir, 'GraphComparison/fig4/fig4d/'];
load([fig4dDir, 'fce.mat']);
fce4d=fce;
load([fig4dDir, 't.mat']);
t4d=t;
load([fig4dDir, 'fig4D.mat']);

%%%% Fig04_Plots
figure
subplot(221)
    plot(t4a,[fce4a],fig4A(:,1),fig4A(:,2));
    title('Fig 4a - forces')
    legend('ce - model','ce - graph')
    xlabel('Time [s]')
    ylabel('Force [N]')
subplot(222)
    plot(t4b,[fce4b],fig4B(:,1),fig4B(:,2))
    title('Fig 4b - forces')
    legend('ce - model','ce - graph')
    xlabel('Time [s]')
    ylabel('Force [N]')
subplot(223)
    plot(t4c,[fce4c],fig4C(:,1),fig4C(:,2))
    title('Fig 4c - forces')
    legend('ce - model','ce - graph')
    xlabel('Time [s]')
    ylabel('Force [N]')
subplot(224)
    plot(t4d,[fce4d],fig4D(:,1),fig4D(:,2))
    title('Fig 4d - forces')
    legend('ce - model','ce - graph')
    xlabel('Time [s]')
    ylabel('Force [N]')

%% Fig05
%%%% Fig05_A_Data
fig5aDir = [matDir, 'GraphComparison/fig5/fig5a/'];
load([fig5aDir, 'fse.mat']);
fse5a=fse;
load([fig5aDir, 'lse.mat']);
lse5a=lse;
load([fig5aDir, 'fig5A.mat']);

%%%% Fig05_B_Data
fig5bDir = [matDir, 'GraphComparison/fig5/fig5b/'];
load([fig5bDir, 'fpe_loop_rel.mat']);
fpe5b=fpe_loop_rel;
load([fig5bDir, 'lpe_length.mat']);
lce5b=lpe_length;
load([fig5bDir, 'fig5B.mat']);

%%%% Fig05_C_Data
fig5cDir = [matDir, 'GraphComparison/fig5/fig5c/'];
load([fig5cDir, 'fce.mat']);
fce5c=f_real;
load([fig5cDir, 'lce.mat']);
lce5c=lce_length;
load([fig5cDir, 'fig5C.mat']);

%%%% Fig05_D_Data
fig5dDir = [matDir, 'GraphComparison/fig5/fig5d/'];
load([fig5dDir, 'fce.mat']);
fce5d=fce;
load([fig5dDir, 'lced.mat']);
lced5d=lced;
load([fig5dDir, 'fig5D.mat']);

%%%% Fig05_Plots
figure;
subplot(221)
    plot(lse5a,fse5a,fig5A(:,1),fig5A(:,2))
    title('Fig 5a')
    xlabel('l_{se} [mm]')
    ylabel('Force SE [N]')
    legend('se - model','se - graph')
    axis([0.0186 0.0197 0 1.3])
subplot(222)
    plot(lce5b,fpe5b,fig5B(:,1)*(10^(-3)),fig5B(:,2))
    title('Fig 5b')
    xlabel('l_{ce} [mm]')
    ylabel('Force PE [N]')
    legend('pe - model','pe - graph')
    %axis([0.01 0.019 0 0.0005])
subplot(223)
    plot(lce5c,fce5c,fig5C(:,1)*(10^(-3)),fig5C(:,2))
    title('Fig 5c')
    xlabel('l_{ce} [mm]')
    ylabel('Force CE [N]')
    legend('ce - model','ce - graph')
subplot(224)
    plot(lced5d,fce5d,fig5D(:,1),fig5D(:,2))
    title('Fig 5d')
    xlabel('l_{ce}-dot [mm]')
    ylabel('Force CE [N]')

%% Fig06
%%%% Fig06_A

%%%% Fig06_B

%%%% Fig06_C



%% Fig07
%%%% Fig07_A

%%%% Fig07_B

%%%% Fig07_C

%%%% Fig07_D