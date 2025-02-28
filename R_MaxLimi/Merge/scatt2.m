%% 20221028 
% 2 dots

function scatt2(data1,data2,yname,c1,c2,cc1,cc2,dotsize,fsize)
nargin;
if  nargin==3
    dotsize=40;
    c1=[1 0.5 0.5];
    c2=[109 161 255]./255;
    cc1=[0.2 0.2 0];
    cc2=[0.2 .4 .5];
    fsize=30;
end

scatter(1,data1,dotsize,'MarkerFaceColor',c1,'MarkerEdgeColor',cc1,'MarkerFaceAlpha',1)
hold on
scatter(2,data2,dotsize,'MarkerFaceColor',c2,'MarkerEdgeColor',cc2,'MarkerFaceAlpha',1)

xlim([0 3])
set(gca,'XTick',[0 1 2 3]) 
name={'','Leader','Follower',' '};

set(gca,'XTickLabel',name,'FontSize',fsize,'FontWeight','bold')  							
ylabel(yname,'FontSize',fsize,'FontWeight','bold')
set(gca,'box','off','linewidth',1,'FontSize',8)
end
