%%%%%%%%%%%%%%%%%%%% tab_2  Fills table 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load('../data/fig_4/cortical_morphs.mat')
load('../data/tab_2/alldata.mat') 

Permutation=[16,15,6,7,3,2,5,1,4,8,12,13,9,11,10,14]; % Reorder cells

oNames=cell(16,1);
for name_Ind=1:16
    oNames{name_Ind}=CortCol.Names{Permutation(name_Ind)};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of cells
Trees=CortCol.Trees;
oNums=zeros(16,1);

for Num_ind=1:16
   iCell=Trees{Permutation(Num_ind)};
   ntree=length(iCell);
   oNums(Num_ind)=ntree;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of pairs
oPairs=zeros(16,1);

for Pair_ind=1:16
   iCell=Trees{Permutation(Pair_ind)};
   ntree=length(iCell);
   nPair=0;
   for Pair_tree_ind=1:ntree
       nPair=nPair+Tab2{Permutation(Pair_ind)}{Pair_tree_ind}.NPair;
   end
   oPairs(Pair_ind)=nPair;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of efferent synapses
oNeff=cell(16,1);

for Neff_ind=1:16
   iCell=Trees{Permutation(Neff_ind)};
   ntree=length(iCell);
   neff=[];
   for Neff_tree_ind=1:ntree
       neff=[neff ; Tab2{Permutation(Neff_ind)}{Neff_tree_ind}.EffVec];
   end
   nMean=mean(neff);
   nSE=sqrt(var(neff)/length(neff));
   oNeff{Neff_ind}=strcat(num2str(round(nMean,2)),'±',num2str(round(nSE,2)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Efferent nearest-neighbour ratio
oNNR=cell(16,1);
opVals=cell(16,1);

nnrscats=[];
for axInd=1:size(measSpat,1)
    for denInd=1:size(measSpat,2)
        nex=size(measSpat{axInd,denInd},1);
            for iex=1:nex
               nnrscats=[nnrscats ; measSpat{axInd,denInd}(iex,1) , measSpat{axInd,denInd}(iex,2) , axInd];
            end
    
    end
end

for nnr_ind=1:16
   iF=find(nnrscats(:,3)==Permutation(nnr_ind));
   
   innr=nnrscats(iF,1);
   nnrMean=mean(innr);
   nnrSEM=sqrt(var(innr)/length(innr));
   oNNR{nnr_ind}=strcat(num2str(round(nnrMean,2)),'±',num2str(round(nnrSEM,2)));
   
   ip=nnrscats(iF,2);
   pMean=mean(ip);
   pSEM=sqrt(var(ip)/length(ip));
   opVals{nnr_ind}=strcat(num2str(round(pMean,2)),'±',num2str(round(pSEM,2)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of afferent synapses

oNaff=cell(16,1);
for Naff_ind=1:16
   iCell=Trees{Permutation(Naff_ind)};
   ntree=length(iCell);
   naff=[];
   naffL=[];
   naffH=[];
   for Naff_tree_ind=1:ntree
       naff=[naff ; Tab2{Permutation(Naff_ind)}{Naff_tree_ind}.AffVec];
       if length(Tab2{Permutation(Naff_ind)}{Naff_tree_ind}.AffVec)>0
           naffL=[naffL ; Tab2{Permutation(Naff_ind)}{Naff_tree_ind}.cInts(:,1)];
           naffH=[naffH ; Tab2{Permutation(Naff_ind)}{Naff_tree_ind}.cInts(:,2)];
       end
   end
   nMean=mean(naff);
   nSE=sqrt(var(naff)/length(naff));
   
   nLMean=mean(naffL);
   nLSE=sqrt(var(naffL)/length(naffL));
   
   nHMean=mean(naffH);
   nHSE=sqrt(var(naffH)/length(naffH));
   
   oNaff{Naff_ind}=strcat(num2str(round(nMean,2)),'±',num2str(round(nSE,2)),'(',num2str(round(nLMean,2)),'±',num2str(round(nLSE,2)),',',num2str(round(nHMean,2)),'±',num2str(round(nHSE,2)),')');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of dendritic compartments (branch)
oMbranch=cell(16,1);

for Mbr_ind=1:16
   iCell=Trees{Permutation(Mbr_ind)};
   ntree=length(iCell);
   iMbr=zeros(ntree,1);
   for Mbr_tree_ind=1:ntree
       iMbr(Mbr_tree_ind)=Tab2{Permutation(Mbr_ind)}{Mbr_tree_ind}.MBranch;
   end
   mMbr=mean(iMbr);
   seMbr=sqrt(var(iMbr)/length(iMbr));
   oMbranch{Mbr_ind}=strcat(num2str(round(mMbr,2)),'±',num2str(round(seMbr,2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proportional number of dendritic compartments (branch), means etc
oMbranchProp=cell(16,1);
oNCompbranch=cell(16,1);
oMubranch=cell(16,1);


for Mbrprop_ind=1:16
   iCell=Trees{Permutation(Mbrprop_ind)};
   ntree=length(iCell);
   
   MbrProp=[];
   
   NCompbra=[];
   NCompbraL=[];
   NCompbraH=[];
    
   MuBra=[];
   MuBraL=[];
   MuBraH=[];
   
   for Mbrprop_tree_ind=1:ntree
       iNv=Tab2{Permutation(Mbrprop_ind)}{Mbrprop_tree_ind}.AffVec;
       if length(Tab2{Permutation(Mbrprop_ind)}{Mbrprop_tree_ind}.AffVec)>0
           iNvL=Tab2{Permutation(Mbrprop_ind)}{Mbrprop_tree_ind}.cInts(:,1);
           iNvH=Tab2{Permutation(Mbrprop_ind)}{Mbrprop_tree_ind}.cInts(:,2);
       end
       if length(iNv)>0
           Mp=ceil(Tab2{Permutation(Mbrprop_ind)}{Mbrprop_tree_ind}.MpropVec.*Tab2{Permutation(Mbrprop_ind)}{Mbrprop_tree_ind}.MBranch);
           Is=Mp>0;
           Mp=Mp(Is);
           iNv=iNv(Is);
           iNvL=iNvL(Is);
           iNvH=iNvH(Is);
           MbrProp=[MbrProp ; Mp];
          h=0;
          for i_h=1:Mp
              h=h+1/i_h;
          end
           NCompbra=[NCompbra ; iNv./(Mp.*h)];  
           NCompbraL=[NCompbraL ; iNvL./(Mp.*h)];
           NCompbraH=[NCompbraH ; iNvH./(Mp.*h)];
           
           thisL=length(iNv);
           if thisL>0
               for pInd=1:thisL
                   [mu] = MeanNMk(iNv(pInd),Mp(pInd));
                   if isnan(mu)
                       mu=0;
                   end
                   MuBra=[MuBra ; mu/Mp(pInd)];
                   
                   [muL] = MeanNMk(iNvL(pInd),Mp(pInd));
                   if isnan(muL)
                       muL=0;
                   end
                   MuBraL=[MuBraL ; muL/Mp(pInd)];
                   
                   [muH] = MeanNMk(iNvH(pInd),Mp(pInd));
                   if isnan(muH) && isnan(mu)
                       muH=0;
                   elseif isnan(muH)
                       muH=Mp(pInd);
                   end
                   MuBraH=[MuBraH ; muH/Mp(pInd)];
               end
           end
       end
   end
   MpMean=mean(MbrProp);
   MpSE=sqrt(var(MbrProp)/length(MbrProp));
   oMbranchProp{Mbrprop_ind}=strcat(num2str(round(MpMean,2)),'±',num2str(round(MpSE,2)));
   
   NCompMean=mean(NCompbra);
   NCompSE=sqrt(var(NCompbra)/length(NCompbra));
   NCompMeanL=mean(NCompbraL);
   NCompSEL=sqrt(var(NCompbraL)/length(NCompbraL));
   NCompMeanH=mean(NCompbraH);
   NCompSEH=sqrt(var(NCompbraH)/length(NCompbraH));
   oNCompbranch{Mbrprop_ind}=strcat(num2str(round(NCompMean,2)),'±',num2str(round(NCompSE,2)),'(',num2str(round(NCompMeanL,2)),'±',num2str(round(NCompSEL,2)),',',num2str(round(NCompMeanH,2)),'±',num2str(round(NCompSEH,2)),')');
  
   MuMean=mean(MuBra);
   MuSE=sqrt(var(MuBra)/length(MuBra));
   MuMeanL=mean(MuBraL);
   MuSEL=sqrt(var(MuBraL)/length(MuBraL));
   MuMeanH=mean(MuBraH);
   MuSEH=sqrt(var(MuBraH)/length(MuBraH));
   oMubranch{Mbrprop_ind}=strcat(num2str(round(MuMean,2)),'±',num2str(round(MuSE,2)),'(',num2str(round(MuMeanL,2)),'±',num2str(round(MuSEL,2)),',',num2str(round(MuMeanH,2)),'±',num2str(round(MuSEH,2)),')');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of dendritic compartments (attenuation)
oMAtten=cell(16,1);

for MAtten_ind=1:16
   iCell=Trees{Permutation(MAtten_ind)};
   ntree=length(iCell);
   iMAtten=zeros(ntree,1);
   for MAtten_tree_ind=1:ntree
       iMAtten(MAtten_tree_ind)=Tab2{Permutation(MAtten_ind)}{MAtten_tree_ind}.MAtten;
   end
   mMAtten=mean(iMAtten);
   seMAtten=sqrt(var(iMAtten)/length(iMAtten));
   oMAtten{MAtten_ind}=strcat(num2str(round(mMAtten,2)),'±',num2str(round(seMAtten,2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proportional number of dendritic compartments (attenuation), means etc

oMAttenProp=cell(16,1);
oNCompAtten=cell(16,1);
oMuAtten=cell(16,1);


for MAttenprop_ind=1:16
   iCell=Trees{Permutation(MAttenprop_ind)};
   ntree=length(iCell);
   
   MAttenProp=[];
   
   NCompAtten=[];
   NCompAttenL=[];
   NCompAttenH=[];
   
   MuAtten=[];
   MuAttenL=[];
   MuAttenH=[];
   for Mbrprop_tree_ind=1:ntree
       iNv=Tab2{Permutation(MAttenprop_ind)}{Mbrprop_tree_ind}.AffVec;
       if length(Tab2{Permutation(MAttenprop_ind)}{Mbrprop_tree_ind}.AffVec)>0
           iNvL=Tab2{Permutation(MAttenprop_ind)}{Mbrprop_tree_ind}.cInts(:,1);
           iNvH=Tab2{Permutation(MAttenprop_ind)}{Mbrprop_tree_ind}.cInts(:,2);
       end
       
       
       if length(iNv)>0
           Mp=ceil(Tab2{Permutation(MAttenprop_ind)}{Mbrprop_tree_ind}.MpropVec.*Tab2{Permutation(MAttenprop_ind)}{Mbrprop_tree_ind}.MAtten);
           Is=Mp>0;
           Mp=Mp(Is);
           iNv=iNv(Is);
           iNvL=iNvL(Is);
           iNvH=iNvH(Is);
           
           MAttenProp=[MAttenProp ; Mp];
             h=0;
          for i_h=1:Mp
              h=h+1/i_h;
          end
             NCompAtten=[NCompAtten ; iNv./(Mp.*h)];
             NCompAttenL=[NCompAttenL ; iNvL./(Mp.*h)];
             NCompAttenH=[NCompAttenH ; iNvH./(Mp.*h)];
             
           thisL=length(iNv);
           if thisL>0
               for pInd=1:thisL
                   [mu] = MeanNMk(iNv(pInd),Mp(pInd));
                     if isnan(mu)
                       mu=0;
                   end
                   MuAtten=[MuAtten ; mu/Mp(pInd)];
                   
                  
                   [muL] = MeanNMk(iNvL(pInd),Mp(pInd));
                   if isnan(muL)
                       muL=0;
                   end
                   MuAttenL=[MuAttenL ; muL/Mp(pInd)];
                   
                   [muH] = MeanNMk(iNvH(pInd),Mp(pInd));
                   if isnan(muH)
                       muH=0;
                   end
                   MuAttenH=[MuAttenH ; muH/Mp(pInd)];
                   
               end
           end
       end
   end
   
   MpMean=mean(MAttenProp);
   MpSE=sqrt(var(MAttenProp)/length(MAttenProp));
   oMAttenProp{MAttenprop_ind}=strcat(num2str(round(MpMean,2)),'±',num2str(round(MpSE,2)));
   
   NCompMean=mean(NCompAtten);
   NCompSE=sqrt(var(NCompAtten)/length(NCompAtten));
   NCompMeanL=mean(NCompAttenL);
   NCompSEL=sqrt(var(NCompAttenL)/length(NCompAttenL));
   NCompMeanH=mean(NCompAttenH);
   NCompSEH=sqrt(var(NCompAttenH)/length(NCompAttenH));
   oNCompAtten{MAttenprop_ind}=strcat(num2str(round(NCompMean,2)),'±',num2str(round(NCompSE,2)),'(',num2str(round(NCompMeanL,2)),'±',num2str(round(NCompSEL,2)),',',num2str(round(NCompMeanH,2)),'±',num2str(round(NCompSEH,2)),')');
 
 
   MuMean=mean(MuAtten);
   MuSE=sqrt(var(MuAtten)/length(MuAtten));
   MuMeanL=mean(MuAttenL);
   MuSEL=sqrt(var(MuAttenL)/length(MuAttenL));
   MuMeanH=mean(MuAttenH);
   MuSEH=sqrt(var(MuAttenH)/length(MuAttenH));
   oMuAtten{MAttenprop_ind}=strcat(num2str(round(MuMean,2)),'±',num2str(round(MuSE,2)),'(',num2str(round(MuMeanL,2)),'±',num2str(round(MuSEL,2)),',',num2str(round(MuMeanH,2)),'±',num2str(round(MuSEH,2)),')'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ColNames={'Cell_type','Number','Pairs','Neff','NNR','p_value','Naff','Mbr_total','Mbr_prop','N_over_NCompbra','Mu_over_Mbra','MAtten_total','Atten_prop','N_over_NCompAtten','Mu_over_MAtten'};

T = table(oNames,oNums,oPairs,oNeff,oNNR,opVals,oNaff,oMbranch,oMbranchProp,oNCompbranch,oMubranch,oMAtten,oMAttenProp,oNCompAtten,oMuAtten,'VariableNames',ColNames);

writetable(T,'../panels/tab_2/table_2.txt');
writetable(T,'../panels/tab_2/table_2.csv');

save('../panels/tab_2/table_2','T');