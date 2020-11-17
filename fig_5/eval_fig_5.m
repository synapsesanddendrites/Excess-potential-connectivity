addpath('../functions/trees')
addpath('../functions')
start_trees

load('../data/fig_5/raw_morphs.mat')

% Mouse
nSamp=10000;
mouseGrid=zeros(nSamp,5);
parfor samp_ind=1:nSamp
    [evalObj] = cell_pair_eval_fig5(oMouse);
   iVec=[evalObj.sep evalObj.sVol evalObj.axLen evalObj.denLen evalObj.N];
   mouseGrid(samp_ind,:)=iVec;
   samp_ind
end
save('../data/fig_5/fig_5mouse.mat','mouseGrid') 
%% Human

HumanGrid=zeros(nSamp,5);
parfor samp_ind=1:nSamp
    [evalObj] = cell_pair_eval_fig5(oHuman);
   iVec=[evalObj.sep evalObj.sVol evalObj.axLen evalObj.denLen evalObj.N];
 HumanGrid(samp_ind,:)=iVec;
 samp_ind
end
save('../data/fig_5/fig_5human.mat','HumanGrid')
clear

%% 

eval_fig_S5
