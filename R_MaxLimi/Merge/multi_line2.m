%% 20221028 折线图
% 2 line

function multi_line2(x1,y1,x2,y2,yname,xname,c1,c2,linesize,fsize)

if  nargin==4

    xname='Time';
    yname='';
    linesize=2;
    fsize=8;
    c1=[1 0.5 0.5];
    c2=[109 161 255]./255;



elseif  nargin==5

    xname='Time';
    linesize=2;
    fsize=8;
    c1=[1 0.5 0.5];
    c2=[109 161 255]./255;


end
plot(x1,y1,'Color',c1,'LineWidth',linesize);hold on
plot(x2,y2,'Color',c2,'LineWidth',linesize)

set(gca,'ycolor','k');

xlabel(xname,'FontSize',fsize,'FontWeight','bold')
ylabel(yname,'FontSize',fsize,'FontWeight','bold')
set(gca,'box','off','linewidth',1,'FontSize',8)
end

