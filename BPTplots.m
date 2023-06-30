function [fig,ax]=BPTplots(plumestruct,var)

arguments
    plumestruct struct                  % structure containing model output
    var string=["radius" "w" "melt"]    % vector of strings containing variable names to plot
end

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
h=abs(diff(r(1:2)))-.05;
cc=repmat(c(1:end-1),[1,m]);
rr=repelem(r(2:end),n);
wh=repmat([w h],[m*n,1]);
pos=[cc' rr' wh];

fn=fieldnames(plumestruct); 

fig=figure('Position',[10 10 200*n+300 300*m+300]);
for i=1:length(var)
    fi=find(fn==var(i));
    ax(i)=axes('Position',pos(i,:));
    plot(plumestruct.(fn{fi}),plumestruct.depth,'LineWidth',2)
    hold on
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
        legend([mh(1) nd(1)],'location','southeast')
    end
end
linkaxes(ax,'y')



