clc
clear
% change
path=pwd;

%% CHANGE %%%%%%%%%%%%%%%%%%
% get basic parameter
ti_num=50;
runsfile=dir([path,'/','FAF*']);
load([path,'/',runsfile(1).name])
% runs_num=floor(length(runsfile)/ti_num);

% capture runs_num for on going data (merge what we have)
run_numbers = [];  % Collect all run numbers

for i = 1:length(runsfile)
    str = runsfile(i).name;
    
    % Use regex to extract the number between 'FAF__run' and '_ti'
    tokens = regexp(str, 'FAF__run(\d+)_ti', 'tokens');
    
    if ~isempty(tokens)
        num = str2double(tokens{1}{1});
        run_numbers(end+1) = num;  % Append to list
    end
end

% Find the maximum run number
runs_num = max(run_numbers);

%%
BEEffort_LT=zeros(T_all/skip,length(leader_value),a_mu);
BEEffort_FT=zeros(T_all/skip,length(leader_value),a_mu);
SmeanT=zeros(T_all/skip,length(leader_value),a_mu);      % average group strength
Effort_LT=zeros(T_all/skip,length(leader_value),a_mu);     % averaghe efforts by rank
Effort_FT=zeros(T_all/skip,length(leader_value),a_mu);    % averaghe efforts by rank
LR_LT=zeros(T_all/skip,length(leader_value),a_mu);
LR_FT=zeros(T_all/skip,length(leader_value),a_mu);
Wmean_LT=zeros(T_all/skip,length(leader_value),a_mu);
Wmean_FT=zeros(T_all/skip,length(leader_value),a_mu);
P_averageT=zeros(T_all/skip,length(leader_value),a_mu);
Wmean_L1T=zeros(T_all/skip,length(leader_value),a_mu);
Wmean_F1T=zeros(T_all/skip,length(leader_value),a_mu);

%%
for runs=1:runs_num
    if runs==123
        
    else
        timesfile=dir([path,'/','FAF__run',num2str(runs),'*']);
        %   timesfile=dir([path,'/','FAF_omg_run1_*']);
        
        %% change
        for ttime=1:length(timesfile)
            load([path,'/',timesfile(ttime).name])
            t_range=T_all/ti_num/skip*(ti-1)+1:T_all/ti_num/skip*ti;
            %% t_range=T_all/20/skip*(ti-1)+1:T_all/20/skip*ti;
            BEEffort_LT(t_range,:,:)=BEMEAN_L;
            BEEffort_FT(t_range,:,:)=BEMEAN_F;
            SmeanT(t_range,:,:)=SMEAN;      % average group strength
            Effort_LT(t_range,:,:)=XMEAN_L;     % averaghe efforts by rank
            Effort_FT(t_range,:,:)=XMEAN_F;    % averaghe efforts by rank
            LR_LT(t_range,:,:)=LRMEAN_L;
            LR_FT(t_range,:,:)=LRMEAN_F;
            Wmean_LT(t_range,:,:)=WMEAN_L;
            Wmean_FT(t_range,:,:)=WMEAN_F;
            P_averageT(t_range,:,:)=GROUP_P;
            Wmean_L1T(t_range,:,:)=WMEAN_L1;
            Wmean_F1T(t_range,:,:)=WMEAN_F1;
            %Round_NumT(t_range,:,:)=Round_Num;
        end
        name=sprintf('FAF_R200_Cost30_run%d.a_mu%d.n%d.b%d.G%d.T%d.w0%d.Lv%d',runs,a_mu,n,b,G,T_all,w_base,length(leader_value));
        mat_file=strrep(name,'.','_');
        mkdir([path,'/Merge/'])
        cd([path,'/Merge/'])
        savepath=[mat_file,'.mat'];
        eval(['save ', savepath,'    TLlist T_all ti_num ti a_mu n G T_range b w_base c0 leader_value skip BEEffort_LT BEEffort_FT SmeanT Effort_LT Effort_FT LR_LT LR_FT Wmean_LT Wmean_FT  Wmean_L1T Wmean_F1T P_averageT'])
        cd(path)
        
        fprintf(['[merging ',num2str(length(timesfile)),'ti files in run ',num2str(runs),']'])
        
    end
end