%% 2023/12/29 Lingyu
% for 2 gene
% There should only one file in result
% if 1,3 plot dotted line too hard to distinguish, change distinguish
% legend update

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

EINL(:,:) = squeeze(mean(mean(INVEST(end,1,:,:,:,:), 6),3)); % 200 x k_level
% Evolution result of leader investment during each round, mean of runs, group
EINF(:,:) = squeeze(mean(mean(mean(INVEST(end,2:end,:,:,:,:), 6),3),2));
% 200 x k_level
% Evolution result of follower investment during each round, mean of runs, group

k5  = [224,136,128]./255;
k20  = [241,188,187]./255;
k200  = [223,225,121]./255;

highc=[0 0 0];
lowc=[0.7 0.7 0.7];

%%
subplot(2,3,1)
l1 = plot(x,XL(:,1),'Color',k5,'LineWidth',4,'LineStyle','-','DisplayName','Leader in k=5');
hold on
l2 = plot(x,XF(:,1),'Color',k5,'LineWidth',4,'LineStyle',':','DisplayName','Follower in k=5');
hold on
l3 = plot(x,XL(:,2),'Color',k20,'LineWidth',4,'LineStyle','-','DisplayName','Leader in k=20');
hold on
l4 = plot(x,XF(:,2),'Color',k20,'LineWidth',4,'LineStyle',':','DisplayName','Follower in k=20');
hold on
l5 = plot(x,XL(:,3),'Color',k200,'LineWidth',4,'LineStyle','-','DisplayName','Leader in k=200');
hold on
l6 = plot(x,XF(:,3),'Color',k200,'LineWidth',4,'LineStyle',':','DisplayName','Follower in k=200');
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

plot(1:5, EINL(1:5,1),'LineWidth',2,'Color',k5,'LineStyle','-')
hold on
plot(1:5, EINF(1:5,1),'LineWidth',2,'Color',k5,'LineStyle',':')
hold on
ylabel('Investment','FontSize',16)
xlabel('Round','FontSize',16)
set(gca,'FontSize',16)
hold off

subplot(2,3,5)

plot(1:20, EINL(1:20,2), 'LineWidth', 2,'Color',k20,'LineStyle','-')
hold on
plot(1:20, EINF(1:20,2), 'LineWidth', 2,'Color',k20,'LineStyle',':')
hold on
ylabel('Investment','FontSize',16)
xlabel('Round','FontSize',16)
set(gca,'FontSize',16)
hold off

subplot(2,3,6)
plot(1:200, EINL(1:200,3), 'LineWidth', 2,'Color',k200,'LineStyle','-')
hold on
plot(1:200, EINF(1:200,3), 'LineWidth', 2,'Color',k200,'LineStyle',':')
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

lgd = legend([l1,l2,l3,l4,l5,l6]);
lgd.Position = [-0.44,0.02,1,0.2];
lgd.Box = 'off';