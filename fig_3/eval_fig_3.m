addpath('../functions/trees')
addpath('../functions')
start_trees


load('../data/fig_3/raw_morphs.mat')
nSamp=10000;

SynPairs=zeros(nSamp,2);
parfor i=1:nSamp
    [output]=cell_pair_eval_fig3(axons,dendrites,2.5)
    SynPairs(i,:)=[output.nEst output.nSyn];
    i
end
SynPairs(isnan(SynPairs))=0;

save('../data/fig_3/alldata.mat','SynPairs') 
clear

%% 

addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/fig_3/raw_morphs.mat')

%------------ Lengths and volumes ------------------------
nDend=length(dendrites);

DendriteLV=zeros(nDend,2);
for dend_ind=1:nDend
    intree=dendrites{dend_ind};
    
    ilen=len_tree(intree);
    iL=sum(ilen(:));
    
    [bound] = boundary_tree(intree);
    iV=bound.V;
    
    DendriteLV(dend_ind,:)=[iL , iV];
end

nAx=length(axons);

AxonLV=zeros(nAx,2);
for ax_ind=1:nAx
    intree=axons{ax_ind};
    
    ilen=len_tree(intree);
    iL=sum(ilen(:));
    
    [bound] = boundary_tree(intree);
    iV=bound.V;
    
   AxonLV(ax_ind,:)=[iL , iV];
end

save('../data/fig_3/supdata.mat','AxonLV','DendriteLV') 
clear


%% Negative parameters

addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/fig_3/alldata.mat')

[PredGrid,MeanGrid,seGrid,sdGrid,sdseGrid,pCons]=collate_synEst(SynPairs);
paramVec=zeros(length(PredGrid),3);
allVec=cell(length(PredGrid),1);
parfor N_ind=2:length(PredGrid)
    N=PredGrid(N_ind)
    [outParams,old_vals]=get_neg_Params(N);
    paramVec(N_ind,:)=[outParams.D outParams.K outParams.r];
    allVec{N_ind}=old_vals;
end
paramVec=smooth_params(allVec);

extraps(1,:) = neg_hyp_tail_fit(PredGrid(PredGrid>3), paramVec(PredGrid>3,1));
extraps(2,:) = neg_hyp_tail_fit(PredGrid(PredGrid>3), paramVec(PredGrid>3,2));
extraps(3,:) = neg_hyp_tail_fit(PredGrid(PredGrid>3), paramVec(PredGrid>3,3));

save('../data/fig_3/confdata.mat','paramVec','extraps') 
clear


%% ; % Confidence intervals;
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/fig_3/confdata.mat','paramVec')

CVec=[25,50,75,95];
nres=10000;
Llines=zeros(length(CVec),nres);
Ulines=zeros(length(CVec),nres);

rawGrid=linspace(0,10,nres);
for c_ind=2:nres
     N=rawGrid(c_ind);
    for val_ind=1:length(CVec)
        val=CVec(val_ind);
        [cint]=conf_int(N,val,paramVec);
        
        Llines(val_ind,c_ind)=cint(1);
        Ulines(val_ind,c_ind)=cint(2);      
    end
end
for c_ind=2:nres
    for val_ind=1:length(CVec)
        if Llines(val_ind,c_ind)<Llines(val_ind,c_ind-1)
            Llines(val_ind,c_ind)=Llines(val_ind,c_ind-1);
        end
        if Ulines(val_ind,c_ind)<Ulines(val_ind,c_ind-1)
            Ulines(val_ind,c_ind)=Ulines(val_ind,c_ind-1);
        end   
    end
end
save('../data/fig_3/conf_int.mat','CVec','rawGrid','Llines','Ulines')
clear
load('../data/fig_3/confdata.mat','paramVec')

CVec=[25,50,75,95];
nres=10000;
Llines=zeros(length(CVec),nres);
Ulines=zeros(length(CVec),nres);

rawGrid=linspace(0,500,nres);
for c_ind=2:nres
     N=rawGrid(c_ind);
    for val_ind=1:length(CVec)
        val=CVec(val_ind);
        [cint]=conf_int(N,val,paramVec);
        
        Llines(val_ind,c_ind)=cint(1);
        Ulines(val_ind,c_ind)=cint(2);      
    end
end
for c_ind=2:nres
    for val_ind=1:length(CVec)
        if Llines(val_ind,c_ind)<Llines(val_ind,c_ind-1)
            Llines(val_ind,c_ind)=Llines(val_ind,c_ind-1);
        end
        if Ulines(val_ind,c_ind)<Ulines(val_ind,c_ind-1)
            Ulines(val_ind,c_ind)=Ulines(val_ind,c_ind-1);
        end   
    end
end
save('../data/fig_3/Lconf_int.mat','CVec','rawGrid','Llines','Ulines')



%% KL Divergence panel 
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/fig_3/confdata.mat','paramVec','extraps')

Ns=10.^linspace(-1,3,1000);
KLdivs=zeros(length(Ns),1);
for N_ind=1:length(Ns)
    KLdivs(N_ind) = KL_Div(Ns(N_ind),paramVec,extraps);
end
save('../data/fig_3/KLdata.mat','Ns','KLdivs')
clear