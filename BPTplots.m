function [fig,ax]=BPTplots(exp,var,amb)

arguments
    exp struct                          % structure containing model output
    var string=["radius" "w" "melt"]    % vector of strings containing variable names to plot
    amb logical=1                       % select whether to include ambient ocean profile in appropriate plots
end

% BPTPLOTS plots outputs of BPTmodel
%   [fig,ax] = BPTplots(exp) From output structure exp from BPTmodel, creates a figure with subplots of the three default 
%   variables (radius, w, and melt). Depths of neutral density and maximum height will also be marked.
%
%   [fig,ax] = BPTplots(exp,var) Creates a figure with subplots of variables listed in var, a string vector of desired 
%   variables. Variables may be any number of the following: radius, w, temp, salt, density, volumeFlux, momentumFlux, and 
%   melt. The term "all" may be substituted to plot all 8. Ambient ocean profile will also be plotted with temp, salinity, 
%   and density.
%   
%   [fig,ax] = BPTplots(exp,var,amb) Use positional argument amb (logical) to plot without ambient profile. 
%
% Examples:
%   % minimum inputs result in plots of radius, w, and melt
%   [fig,ax] = BPTplots(exp)
%
%   % plot all available variables
%   [fig,ax] = BPTplots(exp,'all')
%
%   % plot custom set of variables
%   [fig,ax] = BPTplots(exp,["w", "melt", "temp", "salt"])
%
%   % plot all available variables without ambient profile
%   [fig,ax] = BPTplots(exp,'all',0)
%
% Variables that can be plotted:
%   radius          w               temp            salt        
%   density         volumeFlux      momentumFlux    melt
%
% See also BPTmodel



if strcmp(var,"all")
    var=["radius" "w" "temp" "salt" "density" "volumeFlux" "momentumFlux" "melt"];
end

% determine # rows and columns needed for # var input
m=floor(sqrt(length(var)));
n=ceil(length(var)/m);

%set up axes layout
c=linspace(.07+(4-n)*.01,.99,n+1);
r=linspace(1,.08,m+1);
w=diff(c(1:2))-.03;
h=abs(diff(r(1:2)))-.06;
cc=repmat(c(1:end-1),[1,m]);
rr=repelem(r(2:end),n);
wh=repmat([w h],[m*n,1]);
pos=[cc' rr' wh];

plumestruct=exp.plume;
fn=fieldnames(plumestruct); 
ambstruct=exp.ambient_ocean;
an=fieldnames(ambstruct);
ii=[];

fig=figure('Position',[10 10 200*n+300 300*m+300]);
for i=1:length(var)
    fi=find(fn==var(i));
    ax(i)=axes('Position',pos(i,:));
    h=[]; 
    plot(plumestruct.(fn{fi}),plumestruct.depth,'LineWidth',2)
    hold on
    if amb==1
        if strcmp(var(i),"temp") || strcmp(var(i),"salt") || strcmp(var(i),"density")
            ai=find(an==var(i));
            a=plot(ambstruct.(an{ai}),ambstruct.depth,'LineWidth',2,'Color','k','DisplayName','Ambient Profile');
            ii=[ii i];
            if i==1
                h=a;
            end
        end
    end
    nd=yline(plumestruct.neutralDensity,':','DisplayName','Neutral Density');
    mh=yline(plumestruct.maximumH,'r','DisplayName','Maximum Height');
    if rem(i-1,n)==0
        ylabel(plumestruct.units(1))
    else
        yticklabels([])
    end
    xlabel(plumestruct.units(fi))
    ylim([min(plumestruct.depth) max(plumestruct.depth)])
    set(ax(i),'FontSize',16)
    if i==1
        h=[h mh(1) nd(1)];
        legend(h,'location','southeast','Box','off')
    end
    if strcmp(var(i),'momentumFlux') %put exponent in ticklabel so it doesn't interfere with axis label
        ax(i).XAxis.Exponent=0;
        xtickformat(ax(i),'%.0e') 
    end
end
if ~isempty(ii)
    legend(ax(ii(1)),a,'location','southwest','Box','off')
end
linkaxes(ax,'y')



