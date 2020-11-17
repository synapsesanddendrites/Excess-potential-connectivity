%% Example figure
addpath('../functions/trees')
addpath('../functions')
start_trees

clear
load('../data/fig_4/cortical_morphs.mat')

Trees=CortCol.Trees;
nGroups=length(Groups);
spinedis=3;

 treeInds=Groups{52};
nTree=size(treeInds,1);
figure
hold on

for TreeInd=1:nTree
    
    iTree=Trees{treeInds(TreeInd,1)}{treeInds(TreeInd,2)};
    
            
    iR=iTree.rnames;
    for indR=1:length(iR)
        if contains(iR{indR},'Axon')
            axNum=indR;
        end
    end
    axVals=iTree.R==axNum;
    otherVals=iTree.R~=axNum;
    
    iDend=delete_tree (iTree,find(axVals));
    iAxon=delete_tree(iTree,find(otherVals));
  
    iDend.D=iDend.D+1;
    iAxon.D=iAxon.D+1;
    
    hp= plot_tree(iDend, [1 1 1], [0 0 000], [], [], '-b1');
    set(hp,'facecolor',Colours_Cort(treeInds(TreeInd,1),:)/256,'edgecolor','none');
    
    hp= plot_tree(iAxon, [1 1 1], [0 0 000], [], [], '-b1');
    set(hp,'facecolor',Colours_Cort(treeInds(TreeInd,1),:)/512,'edgecolor','none');
end

view(2);
axis equal off tight
hp= line ([50 100], [-100 -100]);
set(hp,'color',[0 0 0],'linewidth',1);

Ns=strcat('Number of trees:',num2str(nTree));
for type_ind=1:nTree
    iTypes=CortCol.Names(treeInds(type_ind,1));
    Ns=strcat(Ns,', ',iTypes);
end
Ns
saveName = '..\panels\fig_4\fig_4a_Example_column';

set(gcf,'renderer','painter','PaperPositionMode','manual');

tprint(saveName,'-SHR -jpg',[16 16]);
set(gcf,'renderer','painter');
tprint(saveName,'-SHR -eps',[16 16]);
tprint(saveName,'-SHR -tif',[16 16]);

 %% Scatter plot
 clear
load('../data/Colours.mat')
load('../data/fig_4/condata.mat')
Permutation=[16,15,6,7,3,2,5,1,4,8,12,13,9,11,10,14];
 
scats=[];

for axInd=1:size(predCon,1)
    for denInd=1:size(predCon,2)
        nex=length(predCon{axInd,denInd});
            for iex=1:nex
               scats=[scats ; measCon{axInd,denInd}(iex) , pi/2*predCon{axInd,denInd}(iex) , denInd , axInd];
            end
    end
end
scats(isnan(scats))=0;

figure
hold on
plot(0:max(scats(:)),0:max(scats(:)),'black','LineWidth',0.5)
scatter(scats(:,1),scats(:,2),4,Colours_Cort(scats(:,3),:)/256,'filled')
xlim([0 , 150]);
ylim([0 , 150]);
set(gca,'ActivePositionProperty','position','XTick',0:50:150,'YTick',0:50:150,'Xticklabel',0:50:150,'Yticklabel',0:50:150,'XMinorTick','on','YMinorTick','on','ticklength',[0.04 0.02],'fontsize',8,'fontname','helvetica','XGrid','off','YGrid','off','tickdir','out');
saveName = strcat('..\panels\fig_4\fig_4b_Scatter');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
  tprint(saveName,'-SHR -tif',[4 4]);

 ShortNames={'SBC-like','eNGC','L23MC','L23NGC','BTC','BPC','DBC','L23BC','ChC','L23Pyr','L5MC','L5NGC','L5BC','HEC','DC','L5Pyr'};
Locs=[(8:-1:1)' , zeros(8,1) ; (8:-1:1)' , 5*ones(8,1)];
figure
hold on
scatter(Locs(:,2),Locs(:,1),6,Colours_Cort(Permutation,:)/256,'filled');
for ind=1:16
    text(Locs(ind,2)+0.35,Locs(ind,1),strcat(ShortNames{ind}),'FontSize',7,'FontName','helvetica')
end
xlim([-1 8])
ylim([0 9])
axis off
 saveName = strcat('..\panels\fig_4\fig_4b_Scatter_Legend');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
 
 %% Cort_col relations by type
clear
load('../data/Colours.mat')
load('../data/fig_4/condata.mat')
Permutation=[16,15,6,7,3,2,5,1,4,8,12,13,9,11,10,14];

cColours=Colours_Cort(Permutation,:)/256;
cRad=8;

PredGrid=zeros(16);
nGrid=zeros(16); % Number of pairs
for axInd=1:size(predCon,1)
    for denInd=1:size(predCon,2)
        nex=length(predCon{Permutation(axInd),Permutation(denInd)});
        for iex=1:nex
            if ~isnan(predCon{Permutation(axInd),Permutation(denInd)}(iex))
                PredGrid(axInd,denInd)=PredGrid(axInd,denInd)+pi/2*predCon{Permutation(axInd),Permutation(denInd)}(iex);
                nGrid(axInd,denInd)=nGrid(axInd,denInd)+1;
            end
        end
    end
end

MeasGrid=zeros(16);
for axInd=1:size(predCon,1)
    for denInd=1:size(predCon,2)
        nex=length(measCon{Permutation(axInd),Permutation(denInd)});
            for iex=1:nex
               MeasGrid(axInd,denInd)=MeasGrid(axInd,denInd)+measCon{Permutation(axInd),Permutation(denInd)}(iex);
            end
    end   
end
PredGrid(find(nGrid))=PredGrid(find(nGrid))./nGrid(find(nGrid));
MeasGrid(find(nGrid))=MeasGrid(find(nGrid))./nGrid(find(nGrid));

maxPred=max(PredGrid(:));
allVals=PredGrid(:);
allVals=sort(unique(allVals));
figure
hold on
for ind_Est=1:16
    for ind_Meas=1:16
        if PredGrid(ind_Est,ind_Meas)>0
            Gsize=interp1(allVals,linspace(0,1,length(allVals)),PredGrid(ind_Est,ind_Meas));
            scatter3(ind_Est,17-ind_Meas,0,32*Gsize,(PredGrid(ind_Est,ind_Meas)),'filled','square')
        end
    end
end

for type_ind=1:16
        scatter([type_ind , -1 ],[ -1 , 17-type_ind ],cRad,cColours(type_ind,:),'filled');
end
set(gca,'ActivePositionProperty','position','XTick',1:16,'YTick',1:16,'Xticklabel',[],'Yticklabel',[],'ticklength',[0.02 0.02],'fontsize',8,'fontname','helvetica','XGrid','on','YGrid','on','tickdir','in','clipping','off');
xlim([0.5 16.5])
ylim([0.5 16.5])
colormap('Jet')
caxis([0 150])
saveName = strcat('..\panels\fig_4\fig_4c_Cortical_pred_types');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);

maxMeas=max(MeasGrid(:));
allVals=MeasGrid(:);
allVals=sort(unique(allVals));
figure
hold on
for ind_Est=1:16
    for ind_Meas=1:16
        if MeasGrid(ind_Est,ind_Meas)>0
            Gsize=interp1(allVals,linspace(0,1,length(allVals)),MeasGrid(ind_Est,ind_Meas));
            scatter3(ind_Est,17-ind_Meas,0,32*Gsize,(MeasGrid(ind_Est,ind_Meas)),'filled','square')
        end
    end
end

for type_ind=1:16
        scatter([type_ind , -1 ],[ -1 , 17-type_ind ],cRad,cColours(type_ind,:),'filled');
end
set(gca,'ActivePositionProperty','position','XTick',1:16,'YTick',1:16,'Xticklabel',[],'Yticklabel',[],'ticklength',[0.02 0.02],'fontsize',8,'fontname','helvetica','XGrid','on','YGrid','on','tickdir','in','clipping','off');
xlim([0.5 16.5])
ylim([0.5 16.5])
colormap('Jet')
caxis([0 150])
saveName = strcat('..\panels\fig_4\fig_4c_Cortical_meas_types');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
  tprint(saveName,'-SHR -tif',[4 4]);
 
cbh=colorbar('v');
set(cbh,'YTick',0:50:150,'fontsize',8,'fontname','helvetica')

saveName = strcat('..\panels\fig_4\fig_4c_Colourbar');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

 %% Nearest- and all-neighbour ratios
 clear
load('../data/Colours.mat')
load('../data/fig_4/spatdata.mat')
                 
nnrscats=[];
for axInd=1:size(measSpat,1)
    for denInd=1:size(measSpat,2)
        nex=size(measSpat{axInd,denInd},1);
            for iex=1:nex
               nnrscats=[nnrscats ; measSpat{axInd,denInd}(iex,1) , measSpat{axInd,denInd}(iex,2) , axInd];
            end
  
    end
end

anrscats=[];
for axInd=1:size(measSpat,1)
    for denInd=1:size(measSpat,2)
        nex=size(measSpat{axInd,denInd},1);
            for iex=1:nex
               anrscats=[anrscats ; measSpat{axInd,denInd}(iex,3) , measSpat{axInd,denInd}(iex,4) , axInd];
            end
    
    end
end

nnrscats(isnan(nnrscats))=0;
anrscats(isnan(anrscats))=0;

figure
hold on
plot(ones(100,1),linspace(0.000001,1,100),'black','LineWidth',0.5)
plot(linspace(0 ,10,100),0.05*ones(100,1),'black','LineWidth',0.5)
scatter(nnrscats(:,1),nnrscats(:,2),6,Colours_Cort(nnrscats(:,3),:)/256,'filled')
set(gca,'ActivePositionProperty','position','XScale', 'log','YScale', 'log','XTick',10.^linspace(-1,1,3),'YTick',10.^linspace(-4,0,5),'ticklength',[0.02 0.02],'fontsize',8,'fontname','helvetica','tickdir','out','XGrid','on','YGrid','on');
xlim([0.1,10])
ylim([1e-4,1])
saveName = strcat('..\panels\fig_4\fig_4d_NNR');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);


figure
hold on
plot(ones(100,1),linspace(0.000001,1,100),'black','LineWidth',0.5)
plot(linspace(0 ,10,100),0.05*ones(100,1),'black','LineWidth',0.5)
scatter(anrscats(:,1),anrscats(:,2),6,Colours_Cort(anrscats(:,3),:)/256,'filled')
set(gca,'ActivePositionProperty','position','XScale', 'log','YScale', 'log','XTick',10.^linspace(-1,1,3),'YTick',10.^linspace(-4,0,5),'ticklength',[0.04 0.02],'fontsize',8,'fontname','helvetica','tickdir','out','XGrid','on','YGrid','on');
xlim([0.1,10])
ylim([1e-4,1])
saveName = strcat('..\panels\fig_4\fig_4d_ANR');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);
 
 NNRSiglow=nnz(nnrscats(:,1)<1 & nnrscats(:,2)<0.05)
 NNRSighigh=nnz(nnrscats(:,1)>1 & nnrscats(:,2)<0.05)
 ANRSiglow=nnz(anrscats(:,1)<1 & anrscats(:,2)<0.05)
 ANRSighigh=nnz(anrscats(:,1)>1 & anrscats(:,2)<0.05)
 totalPairs=length(nnrscats(:,1)) 
 %% N_complete
clear
load('../data/Colours.mat')
load('../data/fig_4/compdata.mat')
load('../data/fig_4/cortical_morphs.mat')

figure
hold on
Trees=CortCol.Trees;
nClasses=length(Trees);
for ind_Class=1:nClasses
    nCell=length(Trees{ind_Class});   
    for ind_Cell=1:nCell
        Struct=Tab2{ind_Class}{ind_Cell};
        nAff=nnz(Struct.AffVec);
        if nAff>0
            scatter(ceil(Struct.MAtten*Struct.MpropVec(find(Struct.AffVec))),Struct.AffVec(find(Struct.AffVec)),4,Colours_Cort(ind_Class,:)/256,'filled','o')
            scatter(ceil(Struct.MBranch*Struct.MpropVec(find(Struct.AffVec))),Struct.AffVec(find(Struct.AffVec)),4,Colours_Cort(ind_Class,:)/256,'filled','d')
        end
    end
end
M=0:30;
eN=zeros(length(M),1);
for iM=M
    ieN=iM*harmonic(iM);
    eN(iM+1)=ieN;
end
plot(M,eN,'black','LineWidth',0.5)
xlim([0 , 50]);
ylim([0 , 100]);
set(gca,'ActivePositionProperty','position','XTick',0:25:50,'YTick',0:50:100,'Xticklabel',0:25:50,'Yticklabel',0:50:100,'XMinorTick','on','YMinorTick','on','ticklength',[0.04 0.02],'fontsize',8,'fontname','helvetica','XGrid','off','YGrid','off','tickdir','out');
saveName = strcat('..\panels\fig_4\fig_4e_NComplete');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
 
 
  %% Example distributions
 
 clear
load('../data/Colours.mat')
load('../data/fig_4/compdata.mat')
addpath('../functions')

M=50;
Nvals=[5 , 10 , 25 , 50 , 100];
pDists=cell(length(Nvals),1);
for Dist_ind=1:length(Nvals)
   pDists{Dist_ind}=DistinctCompDist(Nvals(Dist_ind),M);   
end

figure
hold on
for Dist_ind=1:length(Nvals) 
    iDist=pDists{Dist_ind};
    iDist1=iDist(:,1);
    iDist2=iDist(:,2);
    isSig=iDist2>0.0005;
    stem(iDist1(isSig),iDist2(isSig),'filled','Color',Colours(Dist_ind,:)/256,'LineWidth',0.5,'MarkerSize',1)
end
xlim([0 , 51]);
ylim([0 , 1]);
set(gca,'ActivePositionProperty','position','XTick',0:25:50,'YTick',0:0.5:1,'Xticklabel',0:25:50,'Yticklabel',0:0.5:1,'XMinorTick','on','YMinorTick','on','ticklength',[0.04 0.02],'fontsize',8,'fontname','helvetica','XGrid','off','YGrid','off','tickdir','out');
legend('N=5','N=10','N=25','N=50','N=100','bkgd','boxoff','FontSize',6,'FontName','helvetica');
saveName = strcat('..\panels\fig_4\fig_4f_Exampledistribution');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
  tprint(saveName,'-SHR -tif',[4 4]);

  %% Example means
clear
load('../data/Colours.mat')
load('../data/fig_4/compdata.mat')
addpath('../functions')

aCs=[1,2,3,4];

NGrid=0:150;
Mcomp=[10 , 25 , 50 , 100];
mus=zeros(length(NGrid),length(Mcomp));
for Nmu_ind=1:length(NGrid)
    for Mmu_ind=1:length(Mcomp)
             mus(Nmu_ind,Mmu_ind)= MeanNMk(NGrid(Nmu_ind),Mcomp(Mmu_ind));
    end
end

figure
hold on
for Dist_ind=1:length(Mcomp) 
  plot(NGrid,mus(:,Dist_ind),'Color',Colours(aCs(Dist_ind),:)/256,'LineWidth',0.5)
  plot(NGrid,Mcomp(Dist_ind)*ones(length(NGrid),1),'Color',Colours(aCs(Dist_ind),:)/256,'LineWidth',0.5,'LineStyle','--','HandleVisibility','off')
end
xlim([0 , 150]);
ylim([0 , 100]);
set(gca,'ActivePositionProperty','position','XTick',0:50:150,'YTick',0:50:100,'Xticklabel',0:50:150,'Yticklabel',0:50:100,'XMinorTick','on','YMinorTick','on','ticklength',[0.04 0.02],'fontsize',8,'fontname','helvetica','XGrid','off','YGrid','off','tickdir','out');
clear
legend('M=10','M=25','M=50','M=100','boxoff','FontSize',6,'FontName','helvetica','Location','northwest');

saveName = strcat('..\panels\fig_4\fig_4g_Tolias_Mus_Examp');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
  %% Data means
clear
load('../data/Colours.mat')
load('../data/fig_4/compdata.mat')
load('../data/fig_4/cortical_morphs.mat')
addpath('../functions')

figure
hold on
Trees=CortCol.Trees;
nClasses=length(Trees);
for ind_Class=1:nClasses
    nCell=length(Trees{ind_Class});
    for ind_Cell=1:nCell
        Struct=Tab2{ind_Class}{ind_Cell};
        nAff=nnz(Struct.AffVec);
        if nAff>0
            MAtten=ceil(Tab2{ind_Class}{ind_Cell}.MpropVec.*Tab2{ind_Class}{ind_Cell}.MAtten);
            MBra=ceil(Tab2{ind_Class}{ind_Cell}.MpropVec.*Tab2{ind_Class}{ind_Cell}.MBranch);
            
            Is=MAtten>0;
            MAtten=MAtten(Is);
            iNv=Struct.AffVec(Is);
            
            thisL=length(iNv);
            iMus=zeros(thisL,2);
            if thisL>0
                for pInd=1:thisL
                    [mu] = MeanNMk(iNv(pInd),MAtten(pInd));
                    iMus(pInd,:)=[MAtten(pInd) , mu];
                end
            end
            scatter(iMus(:,1),iMus(:,2),4,Colours_Cort(ind_Class,:)/256,'filled','o')
            
            Is=MBra>0;
            MBra=MBra(Is);
            iNv=Struct.AffVec(Is);
            
            thisL=length(iNv);
            iMus=zeros(thisL,2);
            if thisL>0
                for pInd=1:thisL
                    [mu] = MeanNMk(iNv(pInd),MBra(pInd));
                    iMus(pInd,:)=[MBra(pInd) , mu];
                end
            end
             scatter(iMus(:,1),iMus(:,2),4,Colours_Cort(ind_Class,:)/256,'filled','d')
        end
       
    end 
end
plot(linspace(0,50,100),linspace(0,50,100),'black','LineWidth',0.5)
xlim([0 , 50]);
ylim([0 , 50]);
set(gca,'ActivePositionProperty','position','XTick',0:25:50,'YTick',0:25:50,'Xticklabel',0:25:50,'Yticklabel',0:25:50,'XMinorTick','on','YMinorTick','on','ticklength',[0.04 0.02],'fontsize',8,'fontname','helvetica','XGrid','off','YGrid','off','tickdir','out');
saveName = strcat('..\panels\fig_4\fig_4g_Data_Mu');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);