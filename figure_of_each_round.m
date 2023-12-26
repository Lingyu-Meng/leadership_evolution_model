clc
clear
filepath=pwd;
file=dir([filepath,'/result/FAF*']);
file_path = [filepath,'/result/',file.name];
load(file_path)

subplot(1,3,1)

plot(1:5,mean(INVMEN_L(end,1:5,1,:),4),'LineWidth',2)
hold on
%scatter(reshape(repmat([1:5]',1,5),1,:),reshape(squeeze(INVMEN_L(end,1:5,1,:)),1,:))
plot(1:5,mean(INVMEN_F(end,1:5,1,:),4),'LineWidth',2)
hold on
ylabel('Investment','FontSize',16)
xlabel('Round','FontSize',16)
set(gca,'FontSize',16)
hold off

subplot(1,3,2)

plot(1:20, mean(INVMEN_L(end,1:20,2,:),4), 'LineWidth', 2)
hold on
plot(1:20, mean(INVMEN_F(end,1:20,2,:),4), 'LineWidth', 2)
hold on
ylabel('Investment','FontSize',16)
xlabel('Round','FontSize',16)
set(gca,'FontSize',16)
hold off



subplot(1,3,3)
plot(1:200, mean(INVMEN_L(end,1:200,3,:),4), 'LineWidth', 2)
hold on
plot(1:200, mean(INVMEN_F(end,1:200,3,:),4), 'LineWidth', 2)
hold on
ylabel('Investment','FontSize',16)
xlabel('Round','FontSize',16)
set(gca,'FontSize',16)
hold off

figure_FontSize=6;
set(gcf,'Position',[0 0 1500 1500]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1);
set(gcf,'color','w')