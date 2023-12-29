%% 2023/12/29 Lingyu
% for 10 gene
% There should only one file in result

clc
clear
filepath=pwd;
file=dir([filepath,'/result/FAF*']);
file_path = [filepath,'/result/',file.name];
load(file_path)

%% prepare data
XL(:,:) = mean(XMEAN_L, 3); % leader Initial_Investment, mean of runs 
XF(:,:) = mean(XMEAN_F, 3); % follower Initial_Investment, mean of runs 

LRL(:,:) = mean(LRMEAN_L, 3); % leader Learning_Rate, mean of runs
LRF(:,:) = mean(LRMEAN_F, 3); % follower Learning_Rate, mean of runs

BEL(:,:) = mean(BEMEAN_L, 3); % leader Last_Round_Behavior, mean of runs
BEF(:,:) = mean(BEMEAN_F, 3); % follower Last_Round_Behavior, mean of runs


k5  = [224,136,128]./255;
k20  = [241,188,187]./255;
k200  = [223,225,121]./255;

highc=[0 0 0];
lowc=[0.7 0.7 0.7];

%%
subplot(2,3,1)

plot(1:skip:T,XL(:,1),'Color',k5,'LineWidth',4,'LineStyle','-')
hold on
plot(1:skip:T,XF(:,1),'Color',k5,'LineWidth',4,'LineStyle','-.')
hold on
plot(1:skip:T,XL(:,2),'Color',k20,'LineWidth',4,'LineStyle','-')
hold on
plot(1:skip:T,XF(:,2),'Color',k20,'LineWidth',4,'LineStyle','-.')
hold on
plot(1:skip:T,XL(:,3),'Color',k200,'LineWidth',4,'LineStyle','-')
hold on
plot(1:skip:T,XF(:,3),'Color',k200,'LineWidth',4,'LineStyle','-.')
hold on

xlabel('Time','FontSize',24,'FontWeight','bold')
ylabel('Initial contribution (x)','FontSize',30,'FontWeight','bold')
ylim([0,10])
set(gca,'box','off','linewidth',3,'FontSize',20)
hold off

subplot(2,3,2)

plot(1:skip:T,LRL(:,1),'Color',k5,'LineWidth',4,'LineStyle','-')
hold on
plot(1:skip:T,LRF(:,1),'Color',k5,'LineWidth',4,'LineStyle','-.')
hold on
plot(1:skip:T,LRL(:,2),'Color',k20,'LineWidth',4,'LineStyle','-')
hold on
plot(1:skip:T,LRF(:,2),'Color',k20,'LineWidth',4,'LineStyle','-.')
hold on
plot(1:skip:T,LRL(:,3),'Color',k200,'LineWidth',4,'LineStyle','-')
hold on
plot(1:skip:T,LRF(:,3),'Color',k200,'LineWidth',4,'LineStyle','-.')
hold on

xlabel('Time','FontSize',24,'FontWeight','bold')
ylabel('Learning rate (u)','FontSize',30,'FontWeight','bold')
ylim([0,1])
set(gca,'box','off','linewidth',3,'FontSize',20)

subplot(2,3,3)

plot(1:skip:T,BEL(:,1),'Color',k5,'LineWidth',4,'LineStyle','-')
hold on
plot(1:skip:T,BEF(:,1),'Color',k5,'LineWidth',4,'LineStyle','-.')
hold on
plot(1:skip:T,BEL(:,2),'Color',k20,'LineWidth',4,'LineStyle','-')
hold on
plot(1:skip:T,BEF(:,2),'Color',k20,'LineWidth',4,'LineStyle','-.')
hold on
plot(1:skip:T,BEL(:,3),'Color',k200,'LineWidth',4,'LineStyle','-')
hold on
plot(1:skip:T,BEF(:,3),'Color',k200,'LineWidth',4,'LineStyle','-.')

ylim([0,10])
ylabel('Contribution in last round','FontSize',20,'FontWeight','bold')
xlabel('Time','FontSize',24,'FontWeight','bold')
set(gca,'box','off','linewidth',3,'FontSize',20)

subplot(2,3,4)

plot(1:5,mean(INVMEN_L(end,1:5,1,:),4),'LineWidth',2)
hold on
%scatter(reshape(repmat([1:5]',1,5),1,:),reshape(squeeze(INVMEN_L(end,1:5,1,:)),1,:))
plot(1:5,mean(INVMEN_F(end,1:5,1,:),4),'LineWidth',2)
hold on
ylabel('Investment','FontSize',16)
xlabel('Round','FontSize',16)
set(gca,'FontSize',16)
hold off

subplot(2,3,5)

plot(1:20, mean(INVMEN_L(end,1:20,2,:),4), 'LineWidth', 2)
hold on
plot(1:20, mean(INVMEN_F(end,1:20,2,:),4), 'LineWidth', 2)
hold on
ylabel('Investment','FontSize',16)
xlabel('Round','FontSize',16)
set(gca,'FontSize',16)
hold off



subplot(2,3,6)
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