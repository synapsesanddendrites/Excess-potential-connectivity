load('../data/fig_3/supdata.mat')
load('../data/Colours.mat')

sz=6;
figure
hold on
scatter(DendriteLV(1,1),DendriteLV(1,2),sz,Colours(1,:)/256,'filled');
scatter(AxonLV(1,1),AxonLV(1,2),sz,Colours(2,:)/256,'filled');
nCell=length(AxonLV(:,1));
for cell_ind=2:nCell
    scatter(DendriteLV(cell_ind,1),DendriteLV(cell_ind,2),sz,Colours(1,:)/256,'filled','HandleVisibility','off');
    scatter(AxonLV(cell_ind,1),AxonLV(cell_ind,2),sz,Colours(2,:)/256,'filled','HandleVisibility','off');

end
xlim([0 , 20000]);
ylim([0 , 10^8]);
set(gca,'ActivePositionProperty','position','XTick',0:5000:20000,'YTick',[0,5,10]*10^7,'XTickLabel',0:5:20,'YTickLabel',[0,5,10],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
legend('Dendrite','Axon','boxoff','FontSize',6,'FontName','helvetica','Location','northwest');

saveName = strcat('..\panels\fig_3\fig_S3a_LVrats');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);


%% Fano factor
clear
load('../data/Colours.mat')
load('../data/fig_3/alldata.mat')
col1=Colours(1,:)/256;
a=2.944;
amin=2.769;
amax=3.119;
b=-0.124;
bmax=-0.246;
bmin=-0.001;

poiss_col=[153, 0, 0]/256;

[PredGrid,~,~,sdGrid,sdseGrid]=collate_synEst(SynPairs);
poiss_fano=ones(size(PredGrid));

figure
hold on
shadedErrorBar(PredGrid',a+PredGrid.^(b-1),[abs((a+PredGrid.^(b-1))-(amin+PredGrid.^(bmin-1)))' , abs((a+PredGrid.^(b-1))-(amax+PredGrid.^(bmax-1)))'])
errorbar(PredGrid,sdGrid.^2./PredGrid',sdseGrid./PredGrid','o','MarkerSize',1,'MarkerFaceColor',col1,'Color',col1,'LineWidth',0.5,'LineStyle','none','CapSize',1);
plot(PredGrid,poiss_fano,'Color',poiss_col,'LineStyle','-','LineWidth',0.5);
ylim([0 10])
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,5,10],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_3\fig_S3b_fano');
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -ti f',[4 4]);

 
%% Skewness
clear
load('../data/Colours.mat')
load('../data/fig_3/alldata.mat')
load('../data/fig_3/confdata.mat','paramVec')

a=0.1184;
amin=0.08062;
amax=0.1509;
b=0.7615;
bmax=0.6268;
bmin=0.8962;

poiss_col=[153, 0, 0]/256;
pol_col=[65,105,225]/256;



[PredGrid,MeanGrid,seGrid,sdGrid,sdseGrid,pCons,~,meas_skews]=collate_synEst(SynPairs);
 
figure
hold on
shadedErrorBar(PredGrid',exp(-a*PredGrid+b),[abs(exp(-amin*PredGrid+bmin)-exp(-a*PredGrid+b))' , abs(exp(-amax*PredGrid+bmax)-exp(-a*PredGrid+b))'])
scatter(PredGrid,meas_skews,4,'black','filled');

polya_skews=zeros(size(PredGrid));
poiss_skews=zeros(size(PredGrid));

a=2.944;
b=-0.124;
for N_ind=1:length(PredGrid)
    N=PredGrid(N_ind);
    p=1-1/(a+N^(b-1));
    r=N*(1-p)/p;
    polya_skews(N_ind)=(1+p)/sqrt(p*r);
    poiss_skews(N_ind)=1/sqrt(N);
  
end

plot(PredGrid,poiss_skews,'Color',poiss_col,'LineStyle','-','LineWidth',0.5);
plot(PredGrid,polya_skews,'Color',pol_col,'LineStyle','-','LineWidth',0.5);


ylim([0 6])
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,3,6],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_3\fig_S3c_skew');
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
 
%% Examples
clear
load('../data/Colours.mat')
load('../data/fig_3/alldata.mat')
load('../data/fig_3/confdata.mat','paramVec')

synPair=SynPairs;
 
Grid=zeros(11,11);

nSamp=size(synPair,1);


poiss_col=[153, 0, 0]/256;
pol_col=[65,105,225]/256;
neg_col=[1, 66, 37]/256;


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


 X=0:10;
 X2=linspace(0,10,50);
figure
hold on
N=1;
scatter(0:10,Grid(:,N+1),4,'black','filled');

ppoiss=poiss_Fits(X,Grid(:,N+1));
ppol=polya_Fits(X,Grid(:,N+1));
pneg=neg_Fits(X,Grid(:,N+1),1);
x=linspace(0,10,500);

plot(x,ppoiss,'Color',poiss_col,'LineStyle','-','LineWidth',0.5);
plot(x,ppol,'Color',pol_col,'LineStyle','-','LineWidth',0.5);
plot(x,pneg,'Color',neg_col,'LineStyle','-','LineWidth',0.5);
 ylim([0 0.3])
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',0:0.1:0.3,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_3\fig_S3e_exp1');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);
 
figure
hold on
N=5;
scatter(0:10,Grid(:,N+1),4,'black','filled');

ppoiss=poiss_Fits(X,Grid(:,N+1));
ppol=polya_Fits(X,Grid(:,N+1));
pneg=neg_Fits(X,Grid(:,N+1),2);
x=linspace(0,10,500);

plot(x,ppoiss,'Color',poiss_col,'LineStyle','-','LineWidth',0.5);
plot(x,ppol,'Color',pol_col,'LineStyle','-','LineWidth',0.5);
plot(x,pneg,'Color',neg_col,'LineStyle','-','LineWidth',0.5);
ylim([0 0.2])

set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',0:0.1:0.2,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_3\fig_S3f_exp2');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);
%% Parameters
clear
load('../data/Colours.mat')
load('../data/fig_3/confdata.mat','paramVec')

 
neg_col=[1, 66, 37]/256;
figure
for p_ind=1:3
   subplot(3,1,p_ind)
   plot(linspace(0,10,50),smooth(paramVec(:,p_ind)),'Color',neg_col,'LineWidth',0.5)
   if p_ind<3
       ylim([0 20])
        set(gca,'ActivePositionProperty','position','XTick',[],'YTick',[0,20],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
   else
       ylim([0 3])
       set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,3],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
   end
   box off
end

set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,3,6],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_3\fig_S3g_params');
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
 
%% KL Divergence
clear
load('../data/fig_3/KLdata.mat','Ns','KLdivs')

figure
semilogx(Ns,smooth(KLdivs),'Color','black','LineWidth',0.5)
ylim([0,2])
xlim([0,1000])
box off
set(gca,'ActivePositionProperty','position','YTick',[0,1,2],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_3\fig_S3g_KL');
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
 
%% Long confidence intervals
clear
load('../data/Colours.mat')
load('../data/fig_3/Lconf_int.mat')
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


xlim([0 500])
ylim([0 500])
set(gca,'ActivePositionProperty','position','XTick',[0,250,500],'YTick',[0,250,500],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_3\fig_S3h_Lconf_int');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
  
%% Full dist
clear
load('../data/fig_3/confdata.mat','paramVec')

PredGrid=linspace(0,10,50);

Grid=zeros(11,11);
for N_ind=2:11
    N=N_ind-1;
    
    D=interp1(PredGrid,paramVec(:,1),N);
    K=interp1(PredGrid,paramVec(:,2),N);
    rho=interp1(PredGrid,paramVec(:,3),N);
    
    [pn]=neg_hyp_dist(0:10,D,K,rho);
    
    Grid(N_ind,:)=pn;
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

saveName = strcat('..\panels\fig_3\fig_S3i_Heat');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

cbh=colorbar('h');
set(cbh,'YTick',0:0.1:0.6,'fontsize',8,'fontname','helvetica')

saveName = strcat('..\panels\fig_3\fig_S3i_Colorbar');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);
