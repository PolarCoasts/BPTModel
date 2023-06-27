% 
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
title(ax1(1),'Example 1','HorizontalAlignment','right');

%% Run model as a point plume, specify entrainment coefficient, and plot just radius and vertical velocity
ex2=BPTmodel(ambD,ambT,ambS,gl,q,type='point',alpha=.08);

[fig2,ax2]=BPTplots(ex2.plume,["radius" "w"]);
title(ax2(1),'Example 2','HorizontalAlignment','right');