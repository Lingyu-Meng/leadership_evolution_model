%% 20221028 折线图
% 2 line

function multi_line3(x1,y1,x2,y2,x3,y3,yname,yname2,xname,c1,c2,c3,linesize,linesize3,fsize)

nargin;
if  nargin==7
    xname='Time';
    yname2='';
    linesize=2;
    linesize3=0.5;
    fsize=8;
    c1=[1 0.5 0.5];
    c2=[109 161 255]./255;
    c3=[0.5 0.5 0.5];
end

plot(x1,y1,'Color',c1,'LineWidth',linesize);hold on
ylabel(yname,'FontSize',fsize,'FontWeight','bold')

[AX,Ha,Hb]=plotyy(x2,y2,x3,y3);
set(Ha,'Color',c2,'LineWidth',linesize);
set(Hb,'Color',c3,'LineWidth',linesize3);
%set(AX(1),'yTick',[0:0.25:1.25])  							% 设置左边Y轴的刻度
%set(AX(2),'yTick',[0:350]) 								% 设置右边Y轴的刻度
d2=get(AX(2),'ylabel');

xlabel(xname,'FontSize',fsize,'FontWeight','bold')
set(AX(1),'Ycolor','k')
set(AX(2),'Ycolor',c3)
set(d2,'string',yname2,'fontsize',fsize,'FontWeight','bold','Color','k');


set(gca,'box','off','linewidth',1,'FontSize',8)


end
