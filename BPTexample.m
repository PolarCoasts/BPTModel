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

[fig1,ax1]=BPTplots(ex1.plume,'all');
title(ax1(1),'Line Plume','Position',[0 1.5],'HorizontalAlignment','left');

%% Run model as a point plume, specify entrainment coefficient, and plot just radius and vertical velocity
ex2=BPTmodel(ambD,ambT,ambS,gl,q,type='point',alpha=.08);

[fig2,ax2]=BPTplots(ex2.plume,["radius" "w"]);
title(ax2(1),'Point Plume','Position',[0 1.5],'HorizontalAlignment','left');

%% Run model for stacked ambient melt plumes and plot default properties
noq=1e-10; %small discharge to initiate plume
width=1; %to get results per 1m of terminus length
uhoriz=0.2; %ambient along-terminus velocity

ex3=BPTmodel(ambD,ambT,ambS,gl,noq,W=width,uh=uhoriz,type='ambient');

[fig3,ax3]=BPTplots(ex3.plume);
title(ax3(1),'Ambient Melt Plumes','Position',[0 1.5],'HorizontalAlignment','left');

%% Run independent melt calculation
% point estimate of melt
melt=melt_calc(-20,uhoriz,4,30);

% profile estimate of melt, specifying ice temp
uprof=linspace(.3,.01,length(ambD));
meltprof=melt_calc(ambD,uprof,ambT,ambS,iceT=-10);
% convert to daily melt rate
meltprof=meltprof*24*60*60;

% plot melt profile
figm=figure('Position',[10 10 500 600]);
axm=axes('Position',[.14 .1 .8 .85]);
plot(meltprof,ambD,'LineWidth',2)
ylabel('Depth (m)')
xlabel('Melt Rate (m/day)')
set(axm,'FontSize',16)
title(axm,'Ambient Melt - No Plume')



