%%%%%%%%%%%%%%%%%%%%% Panel A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example morphologies

addpath('../functions/trees')
addpath('../functions')

start_trees
load('../data/Colours.mat')
load('../data/fig_3/raw_morphs.mat')

col1=Colours(1,:)/256;
col2=Colours(2,:)/256;

tree1=dendrites{42};
tree2=axons{23};

tree2=tran_tree (tree2, [75 0 0]);

d1=tree1.D;
d2=tree2.D;

d1M=mean(d1(:));
d2M=mean(d2(:));

d1(d1>(3*d1M))=d1M;
d2(d2>(3*d1M))=d2M;

tree1.D=d1;
tree2.D=d2+1;

tree1=soma_tree(tree1,20);
tree2=soma_tree(tree2,20);

figure
hp= plot_tree(tree1, [1 1 1], [0 0 000], [], [], '-b1');
set(hp,'facecolor',col1,'edgecolor','none');
view(2);
hp= line ([50 100], [-100 -100]);
set(hp,'color',[0 0 0],'linewidth',1);

hp= plot_tree(tree2, [1 1 1], [0 0 000], [], [], '-b1');
set(hp,'facecolor',col2,'edgecolor','none');
view(2);
hp= line ([50 100], [-100 -100]);
set(hp,'color',[0 0 0],'linewidth',1);

axis equal off tight

saveName = '..\panels\fig_3\fig_3a_Example';

set(gcf,'renderer','painter','PaperPositionMode','manual');

tprint(saveName,'-SHR -jpg',[7 12]);
set(gcf,'renderer','painter');
tprint(saveName,'-SHR -eps',[7 12]);
set(gcf,'renderer','opengl');
tprint(saveName,'-SHR -tif',[7 12]);


%% Mean, SD, and connection probability
clear
load('../data/Colours.mat')
load('../data/fig_3/alldata.mat')
load('../data/fig_3/confdata.mat','paramVec')
col1=Colours(1,:)/256;

[PredGrid,MeanGrid,seGrid,sdGrid,sdseGrid,pCons]=collate_synEst(SynPairs);

poiss_col=[153, 0, 0]/256;
pol_col=[65,105,225]/256;
neg_col=[1, 66, 37]/256;

% Mean
figure
hold on
plot(linspace(0,10),linspace(0,10),'black','LineWidth',0.5)
errorbar(PredGrid,MeanGrid,seGrid,'o','MarkerSize',1,'MarkerFaceColor',col1,'Color',col1,'LineWidth',0.5,'LineStyle','none','CapSize',1);
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,5,10],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_3\fig_3b_MeanN');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);

% SD
a=2.944;
amin=2.769;
amax=3.119;
b=-0.124;
bmax=-0.246;
bmin=-0.001;
figure
hold on

fitvar=[0,a*PredGrid(2:end)+PredGrid(2:end).^(b)];
lovar=[0,amin*PredGrid(2:end)+PredGrid(2:end).^(bmin)];
hivar=[0,amax*PredGrid(2:end)+PredGrid(2:end).^(bmax)];

shadedErrorBar(PredGrid',fitvar,[fitvar'-lovar',hivar'-fitvar'])
errorbar(PredGrid,sdGrid.^2,sdseGrid,'o','MarkerSize',1,'MarkerFaceColor',col1,'Color',col1,'LineWidth',0.5,'LineStyle','none','CapSize',1);
plot(linspace(0,10),linspace(0,10),'Color',poiss_col,'LineStyle','-','LineWidth',0.5)
plot(PredGrid',PredGrid.*(1+lambertw(PredGrid.*(PredGrid.^0.5112-1))),'Color',pol_col,'LineStyle','--')

ylim([0 20])
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,10,20],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_3\fig_3c_sdN');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);

 poiss_cons=zeros(size(PredGrid));
 polya_cons=zeros(size(PredGrid));
 neg_cons=zeros(size(PredGrid));
 
for N_ind=1:length(PredGrid)
   N=PredGrid(N_ind);
   poiss_cons(N_ind)=1-exp(-N);
   
   p=1-1/(a+N^(b-1));
   r=N*(1-p)/p;
   polya_cons(N_ind)=1-(1-p)^r;
   
   D=paramVec(N_ind,1);
   K=paramVec(N_ind,2);
   r=paramVec(N_ind,3);
   
   inSk=neg_hyp_dist(0:1,D,K,r);
   neg_cons(N_ind)=1-inSk(1);  
end

neg_cons(1)=0;
neg_cons=smooth(neg_cons);


figure
hold on
shadedErrorBar(PredGrid',1-exp(-PredGrid.^0.5437)',[abs(exp(-PredGrid.^0.5437)'-exp(-PredGrid.^0.5069)') , abs(exp(-PredGrid.^0.5437)'-exp(-PredGrid.^0.5805)')]);%exp(-0.6845*PredGrid)');%,'r-o','markerfacecolor',[0.5,0.5,0.5])
scatter(PredGrid,pCons,6,'MarkerEdgeColor',col1,'LineWidth',0.5);
plot(PredGrid,poiss_cons,'Color',poiss_col,'LineStyle','-','LineWidth',0.5);
plot(PredGrid,polya_cons,'Color',pol_col,'LineStyle','-','LineWidth',0.5);
plot(PredGrid,neg_cons,'Color',neg_col,'LineStyle','-','LineWidth',0.5);
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,0.5,1],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_3\fig_3d_pCon');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);

 %% Heat maps
 clear
load('../data/Colours.mat')
load('../data/fig_3/alldata.mat')
synPair=SynPairs;
 
Grid=zeros(11,11);

nSamp=size(synPair,1);

for ind=1:nSamp
   ind1=round(synPair(ind,1))+1;
   ind2=synPair(ind,2)+1;
   
   if (ind1<=11) && (ind2<=11)
       Grid(ind1,ind2)=Grid(ind1,ind2)+1;
   end
   
end

for ind=1:11
    ind_vec=Grid(:,ind);
    nI=sum(ind_vec(:));
    Grid(:,ind)=Grid(:,ind)/nI;
end
for ind=1:11
    ind_vec=Grid(ind,:);
    nI=sum(ind_vec(:));
end
allVals=Grid(:);
allVals=sort(unique(allVals));
figure
hold on
for ind_Est=1:11
    nEst=ind_Est-1;
    for ind_Meas=1:11
        nMeas=ind_Meas-1;
        if Grid(ind_Est,ind_Meas)>0
            Gsize=interp1(allVals,linspace(0,1,length(allVals)),Grid(ind_Est,ind_Meas));
            scatter3(nEst,nMeas,0,32*Gsize,(Grid(ind_Est,ind_Meas)),'filled','square')
        end
    end
end
colormap('Jet')
xlim([-0.5 10.5])
ylim([-0.5 10.5])
set(gca,'ActivePositionProperty','position','XTick',0:10,'YTick',0:10,'ticklength',[0.02 0.02],'fontsize',8,'fontname','helvetica','XGrid','on','YGrid','on','tickdir','out');

saveName = strcat('..\panels\fig_3\fig_3e_Heat');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);

cbh=colorbar('h');
set(cbh,'YTick',0:0.1:0.6,'fontsize',8,'fontname','helvetica')

saveName = strcat('..\panels\fig_3\fig_3e_Colorbar');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

%% Confidence intervals
clear
load('../data/Colours.mat')
load('../data/fig_3/conf_int.mat')
col1=Colours(1,:)/256;

figure
hold on
x2 = [rawGrid, fliplr(rawGrid)];
alphaVec=[1,0.75,0.5,0.25];

for conf_ind=length(CVec):-1:2

    inBetween1 = [Llines(conf_ind,:), fliplr(Llines(conf_ind-1,:))];
    fill(x2, inBetween1,col1,'FaceAlpha',alphaVec(conf_ind));
    
    inBetween2 = [Ulines(conf_ind-1,:), fliplr(Ulines(conf_ind,:))];
    fill(x2, inBetween2,col1,'FaceAlpha',alphaVec(conf_ind));
    
end
inBetween = [Llines(1,:), fliplr(Ulines(1,:))];
fill(x2, inBetween,col1);
   
for conf_ind=1:length(CVec)    
   plot(rawGrid,Llines(conf_ind,:),'Color',col1); 
   plot(rawGrid,Ulines(conf_ind,:),'Color',col1); 
end
plot(linspace(0,10),linspace(0,10),'black')

xlim([0 10])
ylim([0 15])
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,5,10],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_3\fig_3e_conf_int');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);