% 2023/12/26 Lingyu
% gene restructured, save all gene and investment
% k = 5, 20, 200
% Purpose: gain mean investment of leader and follower in each round in one
% generation as function of k (rounds per generation).
% Model Setting:
% One threat level (rr = 1.25)
% Fight with nonevolutionary rival (controled by threat level)
%
% parfor restructured 
% Inmportant Output:
% INVMEN_L  % T/skip x 200 x runs(5) each mean investment of leader
% INVMEN_F  % T/skip x 200 x runs(5) each mean investment of follower
% when k < 200, remainder will be zero

% changed form FAF_HJ_0122.m

clc
clear all
close all

%% Group Parameters
G = 500;                  % number of groups
n = 10;                   % number of males per group, also number of females per group
% k = 200;                  % rounds number
b = 4;                    % benefit to each group member (if p=1/G and v_i=1/n)
w_base = 15;              % basal fitness

%% Threat related Parameters
a_mu   = 1;               % number of threat level %cannot change in this script
c0     = 1;               % initial cost
cost_e = 0.03;            % cost_coefficient
alpha  = 1;               % cost exponent -cx^alpha or -c(exp(x^alpha)-1)
beta   = 1;               % group strength exponent: S^beta
N      = n*G;

%% mutation parameter
T          = 20000;       % generation
mu         = 0.05;        % genetic mutation probability
sigma_base = 0.01;        % st. dev. of genetic mutation effect
sigma_lr   = 0.005;       % st. dev. of genetic mutation learning rate
sig        = 0.2;         % Coefficient of variation of attacker effect
skip       = 100;         % number of gens to skip for statistics
save_files = 1;

%% choosing parameter
gene_lr_positive   = 1;   % if 0, the learing rate can be the negetive value
global_selection   = 0;   % if 1, there is group extinction, within-group selection and female dispersal
learning_when_lose = 0;   % if 1, only update the strategy when losing, otherwise update each round no matter win or lose
CP                 = 0;   % if 1, updata the cost each round when losing
random_dispersal_of_males = 1;
IADC = 0;                 % if 1, use specific funtion of fitness which is  W=w_base+P.*(b*V-Cost.*X.^alpha)

%% Learning parameter
learning_coeffcient=1;

%% responsibility
leader_v = 1/n;           % if >1/n, the group will be hierachical, leader will take more benefits
val = [leader_v;repmat((1-leader_v)./(n-1),n-1,1)];
V = reshape(repmat(val,1,G),N,1);

%%
k_set = [5, 20, 200];
runs = 5;
rng('shuffle');
seed=rng;
XMEAN_L  = zeros(T/skip, 3, runs); % generation x k x runs
WMEAN_L  = zeros(T/skip, 3, runs);
LRMEAN_L = zeros(T/skip, 3, runs);
WMEAN_L1 = zeros(T/skip, 3, runs);
XMEAN_F  = zeros(T/skip, 3, runs);
WMEAN_F  = zeros(T/skip, 3, runs);
LRMEAN_F = zeros(T/skip, 3, runs);
GROUP_P  = zeros(T/skip, 3, runs);
WMEAN_F1 = zeros(T/skip, 3, runs);
BEMEAN_L = zeros(T/skip, 3, runs);
BEMEAN_F = zeros(T/skip, 3, runs);
BEMEAN   = zeros(T/skip, 3, runs);
SMEAN    = zeros(T/skip, 3, runs);

INVEST = zeros(T/skip, n, G, 200, 3, runs); % use same scale, each k level record 200 rounds
% generation x subjects x group x round x k_level x runs
GENE_POOL = zeros(T/skip, 2, n, G, 2, 2, 3, runs); 
% generation x leardership x subjects x group x effort/learning x k_level x runs

for k_level=1:3
    %%
    k = k_set(k_level);       % set k
        
    %%
    parfor r = 1:runs % runs
        
        rr = 1.25;              % mean of attacker effect. Actually defined the threat level at here!
        a_sig=sig*rr;               % standard deviation of attacker effect
        
        Gene_pool = zeros(T/skip, 2, n, G, 2, 2); % all alived gene in selected generations
        Invest    = zeros(T/skip, n, G, 200);
        
        BEEffort=zeros(1,T/skip);
        BEEffort_L=zeros(1,T/skip);
        BEEffort_F=zeros(1,T/skip);
        Smean=zeros(1,T/skip);      % average group strength
        Xmean=zeros(1,T/skip);      % average group effort
        Effort=zeros(1,T/skip);     % averaghe efforts by rank
        Effort_L=zeros(1,T/skip);   % averaghe efforts by rank
        Effort_F=zeros(1,T/skip);   % averaghe efforts by rank
        Wmean=zeros(n,T/skip);      % average fitness by rank
        LR_L=zeros(1,T/skip);
        LR_F=zeros(1,T/skip);
        Wmean_L=zeros(1,T/skip);
        Wmean_F=zeros(1,T/skip);
        Xmean_L=zeros(1,T/skip);
        Xmean_F=zeros(1,T/skip);
        P_average=zeros(1,T/skip);
        Wmean1=zeros(n,T/skip);
        Wmean_L1=zeros(1,T/skip);
        Wmean_F1=zeros(1,T/skip);
        
        % initial conditions
        genes = rand(2,n,G,2,2);                        % define the size,(genes,subject,group number,contribution or learning rate,sex)
        genes(:,:,:,1,:) = 2*rr*genes(:,:,:,1,:);       % basic contribution   U(0,2*rr)
        if gene_lr_positive
            genes(:,:,:,2,:) = learning_coeffcient.*genes(:,:,:,2,:);
        else
            genes(:,:,:,2,:) = learning_coeffcient.*(1-2.*genes(:,:,:,2,:));      % learning rate (-1~1)
        end
        
        %*
        fprintf('Rounds t  P  Threat BEEffort_L  BEEffort_F Effort_L Effort_F LearningRate_L LearningRateR_F  Fitness_L  Fitness_F  S  \n');
        cost = repmat(c0,n,G);
        
        for t = 1:T
            W = zeros(N,k);
            x_behavior_update = zeros(n,G);
            follower_lr = zeros(n-1,G);
            leader_lr = zeros(1,G);
            x = zeros(n,G);
            for i = 1:n
                if i == 1
                    leader_lr = reshape(genes(1,1,:,2,2),1,G);
                    x(1,:) = genes(1,i,:,1,2);    % Basic effort of male leader
                else
                    follower_lr(i-1,:) = reshape(genes(2,i,:,2,2),1,G);
                    x(i,:) = genes(2,i,:,1,2);    % Basic effort of male follower
                end  
            end
            x_behavior = x;                       % n x G

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% random choose to fight %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fail_times=zeros(1,G);
            p_round=zeros(G,k);
            
            Invmean_everyGen = zeros(n, G, 200);
            
            for round=1:k
                % avoid too much effort behavior ( to keep w>0)
                % keep the effort behavior is bigger than zero
                x_behavior=reshape(x_behavior,1,N);
                max_behavior=(w_base+b*val)./cost;
                max_behavior=reshape(max_behavior,1,N);
                max_location=find(x_behavior>max_behavior);
                zero_location=find(x_behavior<0);
                
                x_behavior(max_location)=max_behavior(max_location);
                x_behavior(zero_location)=0.0001;
                
                x_behavior=reshape(x_behavior,n,G);
                %--record x_behavior
                Invmean_everyGen(:,:,round) = x_behavior;
                %-----------------------------------------------------------------%
                att=n*(rr+a_sig*randn(G,k)); %power of attack group
                att(att<0)=0;
                
                S=sum(x_behavior(:,:)).^beta; % group strength
                p=S./(S+att(:,round)');           %probability of win 1xG vectors
                
                group_failness=p<0.5;
                fail_times=fail_times+group_failness;
                f=zeros(1,G);
                % for all males in the whole population (Nx1 vectors)
                X=reshape(x_behavior,N,1);
                P=reshape(repmat(p,n,1),N,1);
                Cost=reshape(cost,N,1);
                
                if IADC
                    W(:,round)=w_base+P.*(b*V-Cost.*X.^alpha);     % NEED TO RECOMPUTE WITH TAX!!!!s
                else
                    W(:,round)=w_base+b.*V.*P-Cost.*X.^alpha;     % NEED TO RECOMPUTE WITH TAX!!!!
                end
                p_round(:,round)=p;
                %% update
                if CP
                    ff=fail_times;
                else
                    ff=f;
                end
                
                cost_updata=c0.*exp(cost_e.*ff);   % cost update
                cost=repmat(cost_updata,n,1);
                
                %% stratgy revision
                x(x<0)=0.001;
                if learning_when_lose
                    x_behavior_update=x_behavior;
                    x_behavior_update(1,find(group_failness))=x_behavior(2:n,find(group_failness)).*leader_lr(1,find(group_failness))+x_behavior(1,find(group_failness)).*(1-leader_lr(1,find(group_failness)));
                    x_behavior_update(2:n,find(group_failness))=repmat((x_behavior(1,find(group_failness))),n-1,1).*follower_lr(1,find(group_failness))+x_behavior(2:n,find(group_failness)).*(1-follower_lr(1,find(group_failness)));
                else
                    x_behavior_update(1,:)=mean(x_behavior(2:n,:),1).*leader_lr(1,:)+x_behavior(1,:).*(1-leader_lr(1,:));
                    x_behavior_update(2:n,:)=repmat((x_behavior(1,:)),n-1,1).*follower_lr(1:n-1,:)+x_behavior(2:n,:).*(1-follower_lr(1:n-1,:));
                end
                %define new behavioral
                x_behavior=x_behavior_update;
            end
            
            W_individual_average=mean(W,2);
            W_individual_average=reshape(W_individual_average,n,G);
            p_average=mean(p_round,2);
            p_average(p_average==0)=0.00001;
            
            if global_selection                 % MAY NOT WORK NOW
                %             Female_genes=reshape(genes(:,:,:,1),n+2,N);
                %             Male_genes  =reshape(genes(:,:,:,2),n+2,N);
                %             Fathers=randp(W,2*N,1);
                %             if random_mothers
                %                 Mothers=randi(N,2*N,1);                % if females chosen randomly, evolves to x=b/n
                %             end
                %             R=randi(2,n+2,2*N)-1;                   % indices for free recombination
                %             C=Female_genes(:,Mothers).*R...
                %                 +Male_genes(:,Fathers).*(1-R);  % free recombination
                %             genes(:,:,:,1)=reshape(C(:,1:N),n+2,n,G);
                %             genes(:,:,:,2)=reshape(C(:,N+1:end),n+2,n,G);
            else
                % group extinction-colonization
                new_groups=randp(p_average,1,G);
                genes=genes(:,:,new_groups,:,:);          % redefine genes and p
                x=x(:,new_groups);
                p=p(:,new_groups);
                cost=cost(:,new_groups);
                W_individual_average=W_individual_average(:,new_groups);
                W_individual_average1=W_individual_average./sum(sum(W_individual_average));
                x_behavior=x_behavior(:,new_groups);
                
                % genes come from the mother or the father?
                Rec=randi(2,n,2*n,G)-1; % 4 kinds genes,2n-subjects, group number)
                
                %%
                for g=1:G
                    female_genesx=genes(:,:,g,1,1);        % group female genes (n+2)xn
                    male_genesx  =genes(:,:,g,1,2);        % group male genes (n+2)xn
                    
                    female_genesr=genes(:,:,g,2,1);        % group female genes (n+2)xn
                    male_genesr  =genes(:,:,g,2,2);        % group male genes (n+2)xn
                    
                    w=W_individual_average(:,g);
                    w1=W_individual_average1(:,g);
                    
                    w(w<=0)=0.000001;
                    w1(w1<=0)=0.000001;
                    
                    fathers=randp(w,2*n,1);             % fathers on 2n offspring
                    
                    % each female has two children
                    mothers=repmat(1:n,1,2);
                    % random the order to free recombination
                    randIndex = randperm(2*n);
                    mothers=mothers(randIndex);
                    
                    R=Rec(:,:,g);                       % indices for free recombination
                    
                    % free recombination
                    Cx=female_genesx(:,mothers).*R...
                        +male_genesx(:,fathers).*(1-R);
                    
                    % free recombination
                    Cr=female_genesr(:,mothers).*R...
                        +male_genesr(:,fathers).*(1-R);
                    
                    % assign sexes and leader randomly
                    Cx=Cx(:,randperm(2*n));
                    Cr=Cr(:,randperm(2*n));
                    genes(:,:,g,1,:)=reshape(Cx,2,n,1,1,2);            % redefine genes
                    genes(:,:,g,2,:)=reshape(Cr,2,n,1,1,2);
                end
                % random dispersal of females
                females=reshape(genes(:,:,:,:,1),n,N,2);
                females=females(:,randperm(N),:);
                genes(:,:,:,:,1)=reshape(females,2,n,G,2,1);
                
                if random_dispersal_of_males
                    % random dispersal of males
                    males=reshape(genes(:,:,:,:,2),n,N,2);
                    males=males(:,randperm(N),:);
                    genes(:,:,:,:,2)=reshape(males,2,n,G,2,1);
                end
            end
            
            %%
            if rem(t,skip)==0
                
                Invest(t/skip, :, :, :) = Invmean_everyGen; % record only selected generations
                
                Gene_pool(t/skip, :, :, :, :, :) = genes;
                effort=mean(mean(mean(mean(genes(:,:,:,1,:)))));
                effort_leader=mean(reshape(genes(1,:,:,1,:),1,2*n*G));
                effort_follower=mean(mean(reshape(genes(2:n,:,:,1,:),1,(n-1)*n*2*G)));
                
                Beeffort=mean(mean(x_behavior));
                Beeffort_leader=mean(x_behavior(1,:));
                Beeffort_follower=mean(mean(x_behavior(2:n,:)));
                
                BEEffort_L(:,t/skip)=Beeffort_leader;
                BEEffort_F(:,t/skip)=Beeffort_follower;
                
                Effort(:,t/skip)=effort;
                Effort_L(:,t/skip)=effort_leader;
                Effort_F(:,t/skip)=effort_follower;
                
                LR_L(:,t/skip)= mean(reshape(genes(1,:,:,2,:),1,2*n*G));
                LR_F(:,t/skip)= mean(reshape(genes(2:n,:,:,2,:),1,2*n*(n-1)*G));
                
                Wr1=reshape(W_individual_average1,n,G);              % reshaped W
                Wmean1(:,t/skip)=mean(Wr1,2);                               % mean of relative fitnesses
                Wmean_L1(:,t/skip)=Wmean1(1,t/skip);
                Wmean_F1(:,t/skip)=mean(Wmean1(2:n,t/skip));
                
                Wr=reshape(W_individual_average,n,G);                % reshaped W
                Wmean(:,t/skip)=mean(Wr,2);                                % mean of relative fitnesses
                Wmean_L(:,t/skip)=Wmean(1,t/skip);
                Wmean_F(:,t/skip)=mean(Wmean(2:n,t/skip));
                
                Smean(t/skip)=mean(S);
                P_average(:,t/skip)=mean(p_average);
                
                %% printing stuff
                if rem(t,skip)==0
                    fprintf('%d ',k);
                    fprintf('%d ',t);
                    fprintf('%5.4f ',mean(p_average));
                    fprintf('%5.4f ',rr);
                    fprintf('   %3.2f  %3.2f     %3.2f  %3.2f    %3.2f  %3.2f    %3.2f\n',Effort_L(t/skip),Effort_F(t/skip),LR_L(t/skip),LR_F(t/skip),Wmean_L(t/skip),Wmean_F(t/skip),Smean(t/skip));
                end
            end
            
            %% Mutation
            c=randn;
            if c>0
                % effort
                gene_efforts=reshape(genes(:,:,:,1,:),1,G*n*n*2);
                mutants_be=find(rand(n*n*G*2,1)<mu);
                gene_efforts(mutants_be)=gene_efforts(mutants_be)+sigma_base*randn(length(mutants_be),1)';
                genes(:,:,:,1,:)=reshape(gene_efforts,n,n,G,1,2);
            else
                % learning rate
                gene_lr=reshape(genes(:,:,:,2,:),1,G*n*n*2);
                mutants_lr=find(rand(n*n*G*2,1)<mu);
                gene_lr(mutants_lr)=gene_lr(mutants_lr)+sigma_lr*randn(length(mutants_lr),1)';
                genes(:,:,:,2,:)=reshape(gene_lr,n,n,G,1,2);
            end
            
            %% prevent negative
            % Basic effect
            genes_effort=reshape(genes(:,:,:,1,:),n,n,G,2);
            genes_effort(genes_effort<0)=0.001;
            genes(:,:,:,1,:)=reshape(genes_effort,n,n,G,2);
            % learning rate
            genes_lr=genes(:,:,:,2,:);
            if gene_lr_positive
                genes_lr(genes_lr<0)=0;
                genes_lr(genes_lr>1)=1;
                genes(:,:,:,2,:)=genes_lr;
            else
                genes_lr(genes_lr<-1)=-1;
                genes_lr(genes_lr>1)=1;
                genes(:,:,:,2,:)=genes_lr;
            end
        end
        GROUP_P(:, k_level, r)  = P_average;
        XMEAN_L(:, k_level, r)  = Effort_L;
        WMEAN_L(:, k_level, r)  = Wmean_L;
        LRMEAN_L(:, k_level, r) = LR_L;
        BEMEAN_L(:, k_level, r) = BEEffort_L;
        XMEAN_F(:, k_level, r)  = Effort_F;
        WMEAN_F(:, k_level, r)  = Wmean_F;
        LRMEAN_F(:, k_level, r) = LR_F;
        BEMEAN_F(:, k_level, r) = BEEffort_F;
        BEMEAN(:, k_level, r)   = BEEffort;
        WMEAN_F1(:, k_level, r) = Wmean_F1;
        WMEAN_L1(:, k_level, r) = Wmean_L1;
        SMEAN(:, k_level, r)    = Smean;
        
        INVEST(:, :, :, k_level, r) = Invest;
        GENE_POOL(:, :, :, :, :, :, k_level, r) = Gene_pool;
    end
end    
%%
path=pwd;
if save_files==1
    mkdir('result')
    name=sprintf('FAF.run%d.a_mu%d.n%d.b%d.G%d.T%d.w0%d.round5_20_200.Lv%3.2f',runs,a_mu,n,b,G,T,w_base,leader_v);
    mat_file=strrep(name,'.','_');
    cd([path,'/result/'])
    savepath=[mat_file,'.mat'];
    eval(['save ', savepath,' n G T b w_base c0 leader_v skip SMEAN XMEAN_L WMEAN_L1 WMEAN_L LRMEAN_L  BEMEAN_L BEMEAN_F XMEAN_F WMEAN_F1 WMEAN_F LRMEAN_F GROUP_P INVEST GENE_POOL'])
    cd(path)
end
