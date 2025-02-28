%% 20221028 
% 1 dots

function scatt1(data1,yname,c1,dotsize,fsize)
nargin;
if  nargin==2
    dotsize=40;
    c1=[0.5 0.5 0.5];
    cc1=[0.1 0.1 0.1];
    fsize=30;
end

scatter(1,data1,dotsize,'MarkerFaceColor',c1,'MarkerEdgeColor',cc1,'MarkerFaceAlpha',1)

xlim([0 2])
set(gca,'XTick',[0 1 2]) 
name={'','Group',' '};

set(gca,'XTickLabel',name,'FontSize',fsize,'FontWeight','bold')  							
ylabel(yname,'FontSize',fsize,'FontWeight','bold')
set(gca,'box','off','linewidth',1,'FontSize',8)
end
