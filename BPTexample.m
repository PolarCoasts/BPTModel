% BPTEXAMPLE contains examples for running BPTmodel in a few configurations
clear

%% This example uses an arbitrary ambient profile (cold/fresh at surface, warm/salty at depth)
ambD=-2:-2:-300;
ambT=3-1./logspace(.1,1.5,length(ambD));
ambS=28-5./logspace(.1,1.8,length(ambD));

%% Define grounding line depth and discharge flux
gl=-290; % (m)
q=200;   % (m^3/s)

%% Run model with default values and plot all plume properties
ex1=BPTmodel(ambD,ambT,ambS,gl,q);

[fig1,ax1]=BPTplots(ex1,'all');
title(ax1(1),'Line Plume','Position',[0 1.5],'HorizontalAlignment','left');

%% Run model as a point plume, specify entrainment coefficient, and plot just radius and vertical velocity
ex2=BPTmodel(ambD,ambT,ambS,gl,q,type='point',alpha=.08);

[fig2,ax2]=BPTplots(ex2,["radius" "w"]);
title(ax2(1),'Point Plume','Position',[0 1.5],'HorizontalAlignment','left');

%% Run model for stacked ambient melt plumes and plot default properties
noq=1e-10; %small discharge to initiate plume
uhoriz=0.2; %ambient along-terminus velocity

ex3=BPTmodel(ambD,ambT,ambS,gl,noq,uh=uhoriz,type='stacked');

[fig3,ax3]=BPTplots(ex3);
title(ax3(1),'Stacked Ambient Melt Plumes','Position',[0 1.5],'HorizontalAlignment','left');

%% Run independent melt calculation
% point estimate of melt
melt=melt_calc(-20,uhoriz,4,30);

% profile estimate of melt, specifying ice temp
uprof=linspace(.3,.01,length(ambD));
meltprof=melt_calc(ambD,uprof,ambT,ambS,iceT=0);
% convert to daily melt rate
meltprof=meltprof*24*60*60;

% plot melt profile
figm=figure('Position',[10 10 800 600]);
axu=axes('Position',[.12 .1 .4 .85]);
plot(uprof,ambD,'LineWidth',2,'Color','k')
ylabel('Depth (m)')
xlabel('Along-ice Velocity (m/s)')
set(axu,'FontSize',16)
axm=axes('Position',[.55 .1 .4 .85]);
plot(meltprof,ambD,'LineWidth',2)
yticklabels([])
xlabel('Melt Rate (m/day)')
set(axm,'FontSize',16)
mt=title(axm,'Ambient Melt - No Plume');
mt.Position(1)=-.01;


