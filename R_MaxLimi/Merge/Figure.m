clc
clear
filepath=pwd;
file=dir([filepath,'/FAF*']);
load([filepath,'/',file(1).name])
%%

rn=1;                   % No.? Run
x=10*TLlist;               % att=n*TLlist, power of attack group

%%

% for ff=1:length(file)
%     load([filepath,'/',file(ff).name])
%     XL(:,ff)=Effort_LT;
%     XF(:,ff)=Effort_FT;
%
%     LRL(:,ff)=LR_LT;
%     LRF(:,ff)=LR_FT;
%
%     BEL(:,ff)=BEEffort_LT;
%     BEF(:,ff)=BEEffort_FT;
%
%     GP(:,ff)=P_averageT;
%     WL(:,ff)=Wmean_LT;
%     WF(:,ff)=Wmean_FT;
%
%     Scon(:,ff)=SmeanT;
% end

%%
for ff=1:length(file)
    load([filepath,'/',file(ff).name])
    XL(:,ff)=Effort_LT(10:10:end);
    XF(:,ff)=Effort_FT(10:10:end);

    LRL(:,ff)=LR_LT(10:10:end);
    LRF(:,ff)=LR_FT(10:10:end);

    BEL(:,ff)=BEEffort_LT(10:10:end);
    BEF(:,ff)=BEEffort_FT(10:10:end);

    GP(:,ff)=P_averageT(10:10:end);
    WL(:,ff)=Wmean_LT(10:10:end);
    WF(:,ff)=Wmean_FT(10:10:end);

    Scon(:,ff)=SmeanT(10:10:end);
 %   RN(:,ff)=Round_NumT(10:10:end);
end





%%

c1=[1 0.5 0.5];c2=[109 161 255]./255;
highc=[0 0 0];lowc=[0.7 0.7 0.7];
hl_c=[0.65 0 0];hf_c=[0 44 126]./255;

%% figure1 [ an example run rn=1 ]
IC_L=XL(:,rn);  IC_F=XF(:,rn);   LR_L=LRL(:,rn); LR_F=LRF(:,rn);
BE_L=BEL(:,rn); BE_F=BEF(:,rn);  FIT_L=WL(:,rn); FIT_F=WF(:,rn);
SR=GP(:,rn);    TCon=Scon(:,rn); %rounum=RN(:,rn);

x=1:skip*10:T_all; xname='Time';
yname2='the number of rounds';
c1=[1 0.5 0.5];c2=[109 161 255]./255;c3=[0.5 0.5 0.5];

% Figure 1a
figure;set(gcf,'Position',[0 40 900 600],'color','w');
subplot(2,1,1);yname='Initial contribution (x)';% X
multi_line2(x,IC_L,x,IC_F,yname)

subplot(2,1,2);yname='Learning rate (u)';% u
multi_line2(x,LR_L,x,LR_F,yname)

% Figure 1b
figure;set(gcf,'Position',[50 50 900 600],'color','w');
subplot(3,1,1);yname='Contribution in last round';% X
multi_line2(x,BE_L,x,BE_F,yname)

subplot(3,1,2);yname='Individual Fitness';% u
multi_line2(x,FIT_L,x,FIT_F,yname)

subplot(3,1,3);yname='Success Rate';% u
multi_line2(x,SR,x,SR,yname)



%% figure2 [ repeated runs ]

figure;set(gcf,'Position',[0 40 900 600],'color','w');
for rn=1:length(file)
    IC_L=XL(:,rn);  IC_F=XF(:,rn);
    LR_L=LRL(:,rn); LR_F=LRF(:,rn);
   % rounum=RN(:,rn);

    subplot(2,1,1);yname='Initial contribution (x)';% X
    multi_line2(x,IC_L,x,IC_F,yname)

    subplot(2,1,2);yname='Learning rate (u)';% u
    multi_line2(x,LR_L,x,LR_F,yname)
end


figure;set(gcf,'Position',[50 50 900 600],'color','w');
for rn=1:length(file)
    BE_L=BEL(:,rn); BE_F=BEF(:,rn);
    FIT_L=WL(:,rn); FIT_F=WF(:,rn);
    SR=GP(:,rn);    TCon=Scon(:,rn); %rounum=RN(:,rn);

    %subplot(4,1,1);
    subplot(3,1,1);
    yname='Contribution in last round';% X
    multi_line2(x,BE_L,x,BE_F,yname,xname,c1,c2,0.5,8)

    %subplot(4,1,2);
    subplot(3,1,2);
    yname='Individual Fitness';% u
    multi_line2(x,FIT_L,x,FIT_F,yname,xname,c1,c2,0.5,8)

    %subplot(4,1,3);
    subplot(3,1,3);
    yname='Success Rate';% u
    multi_line2(x,SR,x,SR,yname,xname,c1,c2,0.5,8)

    %     subplot(4,1,4);yname='Group Contribution';% u
    %     multi_line2(x,TCon,x,TCon,yname,xname,c1,c2,0.5,8)
end


%% figure3 [ result for envolution ]

x=[1,2]; % leader / follower

figure;set(gcf,'Position',[0 40 1800 200],'color','w');

subplot(1,6,1);yname='Initial contribution (x)';
scatt2(XL(end,:),XF(end,:),yname)
ylim([0 1.2])

subplot(1,6,2);yname='Learning rate (u)';
scatt2(LRL(end,:),LRF(end,:),yname)
ylim([0 0.2])

subplot(1,6,3);yname='Individual Fitness';% X
scatt2(WL(end,:),WF(end,:),yname)

subplot(1,6,4);yname='Behavioral';% u
scatt2(BEL(end,:),BEF(end,:),yname)

subplot(1,6,5);yname='Success Rate';% u
scatt1(GP(end,:),yname)
ylim([0 0.5])

subplot(1,6,6);yname='Group Contribution';% u
scatt1(Scon(end,:),yname)
ylim([0 6])