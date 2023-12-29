%% 2023/12/29 Lingyu
% for 10 gene
% There should only one file in result
% if 1,3 plot dotted line too hard to distinguish, change distinguish

clc
clear
filepath=pwd;
file=dir([filepath,'/result/FAF*']);
file_path = [filepath,'/result/',file.name];
load(file_path)

distinguish = 1; % affect 1,3 only now

%% prepare data
x = 1:skip:T;
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

plot(x,XL(:,1),'Color',k5,'LineWidth',4,'LineStyle','-')
hold on
plot(x,XF(:,1),'Color',k5,'LineWidth',4,'LineStyle',':')
hold on
plot(x,XL(:,2),'Color',k20,'LineWidth',4,'LineStyle','-')
hold on
plot(x,XF(:,2),'Color',k20,'LineWidth',4,'LineStyle',':')
hold on
plot(x,XL(:,3),'Color',k200,'LineWidth',4,'LineStyle','-')
hold on
plot(x,XF(:,3),'Color',k200,'LineWidth',4,'LineStyle',':')
hold on

xlabel('Time','FontSize',24,'FontWeight','bold')
ylabel('Initial contribution (x)','FontSize',30,'FontWeight','bold')
ylim([0,10])
set(gca,'box','off','linewidth',3,'FontSize',20)
hold off

subplot(2,3,2)

plot(x,LRL(:,1),'Color',k5,'LineWidth',2,'LineStyle','-')
hold on
plot(x,LRF(:,1),'Color',k5,'LineWidth',2,'LineStyle',':')
hold on
plot(x,LRL(:,2),'Color',k20,'LineWidth',2,'LineStyle','-')
hold on
plot(x,LRF(:,2),'Color',k20,'LineWidth',2,'LineStyle',':')
hold on
plot(x,LRL(:,3),'Color',k200,'LineWidth',2,'LineStyle','-')
hold on
plot(x,LRF(:,3),'Color',k200,'LineWidth',2,'LineStyle',':')
hold on

xlabel('Time','FontSize',24,'FontWeight','bold')
ylabel('Learning rate (u)','FontSize',30,'FontWeight','bold')
ylim([0,1])
set(gca,'box','off','linewidth',3,'FontSize',20)

subplot(2,3,3)
idx = x;
be_l = BEL;
be_f = BEF;
if distinguish
idx(1:2:end) = [];
be_l(1:2:end,:) = [];
be_f(1:2:end,:) = [];
end

plot(idx,be_l(:,1),'Color',k5,'LineWidth',2,'LineStyle','-')
hold on
plot(idx,be_f(:,1),'Color',k5,'LineWidth',2,'LineStyle',':')
hold on
plot(idx,be_l(:,2),'Color',k20,'LineWidth',2,'LineStyle','-')
hold on
plot(idx,be_f(:,2),'Color',k20,'LineWidth',2,'LineStyle',':')
hold on
plot(idx,be_l(:,3),'Color',k200,'LineWidth',2,'LineStyle','-')
hold on
plot(idx,be_f(:,3),'Color',k200,'LineWidth',2,'LineStyle',':')
ylim([0,10])
ylabel('Contribution in last round','FontSize',20,'FontWeight','bold')
xlabel('Time','FontSize',24,'FontWeight','bold')
set(gca,'box','off','linewidth',3,'FontSize',20)

subplot(2,3,4)

plot(1:5,mean(INVMEN_L(end,1:5,1,:),4),'LineWidth',2,'Color',k5,'LineStyle','-')
hold on
%scatter(reshape(repmat([1:5]',1,5),1,:),reshape(squeeze(INVMEN_L(end,1:5,1,:)),1,:))
plot(1:5,mean(INVMEN_F(end,1:5,1,:),4),'LineWidth',2,'Color',k5,'LineStyle',':')
hold on
ylabel('Investment','FontSize',16)
xlabel('Round','FontSize',16)
set(gca,'FontSize',16)
hold off

subplot(2,3,5)

plot(1:20, mean(INVMEN_L(end,1:20,2,:),4), 'LineWidth', 2,'Color',k20,'LineStyle','-')
hold on
plot(1:20, mean(INVMEN_F(end,1:20,2,:),4), 'LineWidth', 2,'Color',k20,'LineStyle',':')
hold on
ylabel('Investment','FontSize',16)
xlabel('Round','FontSize',16)
set(gca,'FontSize',16)
hold off



subplot(2,3,6)
plot(1:200, mean(INVMEN_L(end,1:200,3,:),4), 'LineWidth', 2,'Color',k200,'LineStyle','-')
hold on
plot(1:200, mean(INVMEN_F(end,1:200,3,:),4), 'LineWidth', 2,'Color',k200,'LineStyle',':')
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