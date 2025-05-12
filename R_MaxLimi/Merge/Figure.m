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
for r = 1:a_mu
    for ff=1:length(file)
        load([filepath,'/',file(ff).name])
        XL(:,ff,r)=Effort_LT(10:10:end,r);
        XF(:,ff,r)=Effort_FT(10:10:end,r);
        
        LRL(:,ff,r)=LR_LT(10:10:end,r);
        LRF(:,ff,r)=LR_FT(10:10:end,r);
        
        BEL(:,ff,r)=BEEffort_LT(10:10:end,r);
        BEF(:,ff,r)=BEEffort_FT(10:10:end,r);
        
        GP(:,ff,r)=P_averageT(10:10:end,r);
        WL(:,ff,r)=Wmean_LT(10:10:end,r);
        WF(:,ff,r)=Wmean_FT(10:10:end,r);
        
        Scon(:,ff,r)=SmeanT(10:10:end,r);
        %   RN(:,ff)=Round_NumT(10:10:end);
    end
end





%%

c1=[1 0.5 0.5];c2=[109 161 255]./255;
highc=[0 0 0];lowc=[0.7 0.7 0.7];
hl_c=[0.65 0 0];hf_c=[0 44 126]./255;

% plot figures by a_mu
for r = 1:a_mu
    %% figure1 [ an example run rn=1 ]
    IC_L=XL(:,rn,r);  IC_F=XF(:,rn,r);   LR_L=LRL(:,rn,r); LR_F=LRF(:,rn,r);
    BE_L=BEL(:,rn,r); BE_F=BEF(:,rn,r);  FIT_L=WL(:,rn,r); FIT_F=WF(:,rn,r);
    SR=GP(:,rn,r);    TCon=Scon(:,rn,r); %rounum=RN(:,rn);
    
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
    
    sgtitle(['Threat Level:',num2str(TLlist(r)),' w: ',num2str(w_base),' b: ',...
    num2str(b),' c0: ', num2str(c0)]); % sgtitle in MATLAB (R2018b or newer)
    
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
    
    sgtitle(['Threat Level:',num2str(TLlist(r)),' w: ',num2str(w_base),' b: ',...
    num2str(b),' c0: ', num2str(c0)]);
    
    %% figure3 [ result for envolution ]
    
    x=[1,2]; % leader / follower
    
    figure;set(gcf,'Position',[0 40 1800 200],'color','w');
    
    subplot(1,6,1);yname='Initial contribution (x)';
    scatt2(XL(end,:,r),XF(end,:,r),yname)
    ylim([0 1.2])
    
    subplot(1,6,2);yname='Learning rate (u)';
    scatt2(LRL(end,:,r),LRF(end,:,r),yname)
    ylim([0 0.2])
    
    subplot(1,6,3);yname='Individual Fitness';% X
    scatt2(WL(end,:,r),WF(end,:,r),yname)
    
    subplot(1,6,4);yname='Behavioral';% u
    scatt2(BEL(end,:,r),BEF(end,:,r),yname)
    
    subplot(1,6,5);yname='Success Rate';% u
    scatt1(GP(end,:,r),yname)
    ylim([0 0.5])
    
    subplot(1,6,6);yname='Group Contribution';% u
    scatt1(Scon(end,:,r),yname)
    ylim([0 6])
    
    sgtitle(['Threat Level:',num2str(TLlist(r)),' w: ',num2str(w_base),' b: ',...
    num2str(b),' c0: ', num2str(c0)]);
end