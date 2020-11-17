% Length and volume of cortical data
clear
load('../data/fig_4/supdata.mat') 

sz=6;

Axscats=[];
for ind_Class=1:size(AxCell,1)
    nTree=length(AxCell{ind_Class});
    for ind_Tree=1:nTree       
        Axscats=[Axscats ; AxCell{ind_Class}(ind_Tree,:) , ind_Class];
    end
end
Axscats(isnan(Axscats))=0;



figure
hold on
nCell=length(Axscats);
pCell=randperm(nCell);
for cell_ind=2:nCell
    scatter(Axscats(pCell(cell_ind),1),Axscats(pCell(cell_ind),2),sz,Colours_Cort(Axscats(pCell(cell_ind),3),:)/256,'filled');
end
xlim([0 , 80000]);
ylim([0 , 5*10^7]);
set(gca,'ActivePositionProperty','position','XTick',0:20000:80000,'YTick',[0,2.5,5]*10^7,'XTickLabel',0:20:80,'YTickLabel',[0,2.5,5],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');


saveName = strcat('..\panels\fig_4\fig_S4b_AxonLvsV');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
  tprint(saveName,'-SHR -tif',[4 4]);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Denscats=[];
for ind_Class=1:size(DendCell,1)
    nTree=length(DendCell{ind_Class});
    for ind_Tree=1:nTree        
        Denscats=[Denscats ; DendCell{ind_Class}(ind_Tree,:) , ind_Class];
    end
end
Denscats(isnan(Denscats))=0;



figure
hold on
nCell=length(Denscats);
pCell=randperm(nCell);
for cell_ind=2:nCell
    scatter(Denscats(pCell(cell_ind),1),Denscats(pCell(cell_ind),2),sz,Colours_Cort(Denscats(pCell(cell_ind),3),:)/256,'filled');
end
xlim([0 , 10000]);
ylim([0 , 10^7]);
set(gca,'ActivePositionProperty','position','XTick',0:5000:10000,'YTick',[0,5,10]*10^6,'XTickLabel',0:5:10,'YTickLabel',[0,0.5,1],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');


saveName = strcat('..\panels\fig_4\fig_S4a_DendLvsV');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ratio of all- and nearest-neighbour ratios
clear
load('../data/fig_4/spatdata.mat')
sz=6;

nnranrscats=[];
for axInd=1:size(measSpat,1)
    for denInd=1:size(measSpat,2)
        nex=size(measSpat{axInd,denInd},1);
            for iex=1:nex
               nnranrscats=[nnranrscats ; measSpat{axInd,denInd}(iex,1) , measSpat{axInd,denInd}(iex,3) , axInd];
            end
    
    end
end

figure
hold on
scatter(nnranrscats(:,1),nnranrscats(:,2),sz,Colours_Cort(nnranrscats(:,3),:)/256,'filled')
 plot(linspace(0 , max(nnranrscats(:,1)),100),linspace(0 , max(nnranrscats(:,1)),100),'black')
xlim([0 , max(nnranrscats(:,1))]);
ylim([0 , max(nnranrscats(:,1))]);
set(gca,'ActivePositionProperty','position','XTick',0:5:10,'YTick',0:5:10,'XTickLabel',0:5:10,'YTickLabel',0:5:10,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');


saveName = strcat('..\panels\fig_4\fig_S4c_ANRvsNNR');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
 
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Volume ratios
clear
load('../data/fig_4/spatdata.mat') 
  
volrat=[];
for axInd=1:size(vol_rats,1)
    for denInd=1:size(vol_rats,2)
        nex=size(vol_rats{axInd,denInd},1);
            for iex=1:nex
               volrat=[volrat ; vol_rats{axInd,denInd}(iex,1) , vol_rats{axInd,denInd}(iex,2) , axInd];
            end
    
    end
end
  

figure
hold on
plot(ones(100,1),linspace(0.000001,1,100),'black','LineWidth',0.5)
plot(linspace(0 ,10,100),0.05*ones(100,1),'black','LineWidth',0.5)
scatter(volrat(:,1),volrat(:,2),4,Colours_Cort(volrat(:,3),:)/256,'filled')
set(gca,'ActivePositionProperty','position','XScale', 'log','YScale', 'log','XTick',10.^linspace(-1,1,3),'YTick',10.^linspace(-4,0,5),'ticklength',[0.02 0.02],'fontsize',8,'fontname','helvetica','tickdir','out','XGrid','on','YGrid','on');
xlim([0.1 10])
ylim([1e-4 1])
saveName = strcat('..\panels\fig_4\fig_S4d_VolRat');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Relationship between estimates of compartment number
clear
sz=4;
load('../data/fig_4/compdata.mat')
load('../data/fig_4/cortical_morphs.mat')
Trees=CortCol.Trees;

MBraVSMAtten=[];
for ind_Class=1:length(Tab2)
    for ind_Cell=1:length(Tab2{ind_Class})
        MBraVSMAtten=[MBraVSMAtten ; Tab2{ind_Class}{ind_Cell}.MBranch , Tab2{ind_Class}{ind_Cell}.MAtten , ind_Class];
    end
end

figure
hold on
for cell_ind=1:length(MBraVSMAtten)
    scatter(MBraVSMAtten(cell_ind,1),MBraVSMAtten(cell_ind,2),sz,'MarkerEdgeColor',Colours_Cort(MBraVSMAtten(cell_ind,3),:)/256,'MarkerFaceColor',Colours_Cort(MBraVSMAtten(cell_ind,3),:)/256,'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
end

oMbranch=zeros(16,2);
for Mbr_ind=1:16
   iCell=Trees{Mbr_ind};
   ntree=length(iCell);
   iMbr=zeros(ntree,1);
   for Mbr_tree_ind=1:ntree
       iMbr(Mbr_tree_ind)=Tab2{Mbr_ind}{Mbr_tree_ind}.MBranch;
   end
   mMbr=mean(iMbr);
   seMbr=sqrt(var(iMbr)/length(iMbr));
   oMbranch(Mbr_ind,1)=mMbr;
   oMbranch(Mbr_ind,2)=seMbr;
end

oMAtten=zeros(16,2);
for MAtten_ind=1:16
   iCell=Trees{MAtten_ind};
   ntree=length(iCell);
   iMAtten=zeros(ntree,1);
   for MAtten_tree_ind=1:ntree
       iMAtten(MAtten_tree_ind)=Tab2{MAtten_ind}{MAtten_tree_ind}.MAtten;
   end
   mMAtten=mean(iMAtten);
   seMAtten=sqrt(var(iMAtten)/length(iMAtten));
   oMAtten(MAtten_ind,1)= mMAtten;
   oMAtten(MAtten_ind,2)=seMAtten; 
end
[~,brI] = sort(oMbranch(:,1));

for type_ind=1:16
    iInd=brI(type_ind);
    errorbar(oMbranch(iInd,1),oMAtten(iInd,1),oMAtten(iInd,2),oMAtten(iInd,2),oMbranch(iInd,2),oMbranch(iInd,2),'o','MarkerSize',0.5*sz,'MarkerEdgeColor',Colours_Cort(iInd,:)/256,'MarkerFaceColor',Colours_Cort(iInd,:)/256,'Color',Colours_Cort(iInd,:)/256,'CapSize',sz)
end
xlim([0 , 80]);
ylim([0 , 50]);
set(gca,'ActivePositionProperty','position','XTick',0:40:80,'YTick',0:25:50,'XTickLabel',0:40:80,'YTickLabel',0:25:50,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');

  saveName = strcat('..\panels\fig_4\fig_S4e_MbraVSMatten');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
%  
%% Number of compartments (attenuation) as a function of threshold
clear
load('../data/fig_4/supdata.mat') 

figure
hold on
for cell_ind=1:16
errorbar(thresVec,AllMeans(cell_ind,:),AllSEMs(cell_ind,:),'o','MarkerSize',1,'MarkerFaceColor',Colours_Cort(cell_ind,:)/256,'Color',Colours_Cort(cell_ind,:)/256,'LineWidth',0.5,'LineStyle','--','CapSize',1);
end
set(gca,'ActivePositionProperty','position','XTick',[0,0.5,1],'YTick',[0,10,20,30],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_4\fig_S4f_MattenbyThres');   
xlim([0,1])
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]); 

%% Legend
clear
load('../data/fig_4/supdata.mat')
Permutation=[16,15,6,7,3,2,5,1,4,8,12,13,9,11,10,14];
 ShortNames={'SBC-like','eNGC','L23MC','L23NGC','BTC','BPC','DBC','L23BC','ChC','L23Pyr','L5MC','L5NGC','L5BC','HEC','DC','L5Pyr'};
Locs=[(4:-1:1)' , zeros(4,1) ; (4:-1:1)' , 7*ones(4,1); (4:-1:1)' , 14*ones(4,1); (4:-1:1)' , 21*ones(4,1)];
figure
hold on
scatter(Locs(:,2),Locs(:,1),6,Colours_Cort(Permutation,:)/256,'filled');
for ind=1:16
    text(Locs(ind,2)+0.35,Locs(ind,1),strcat(ShortNames{ind}),'FontSize',7,'FontName','helvetica')
end
xlim([-1 25])
ylim([0 5])
axis off
 saveName = strcat('..\panels\fig_4\fig_S4Legend');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[8 2]);
 tprint(saveName,'-SHR -eps',[8 2]);
 tprint(saveName,'-SHR -tif',[4 4]);