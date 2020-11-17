%%%%%%%%%%%%%%%%%%%%% Panel A1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shape vec morphologies

clear
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/Colours.mat')
load('../data/fig_2/morphs_fig_2a.mat')

colDilution=[1 , 0.8, 0.6 , 0.4];

sphere_tree.D=sphere_tree.D+1;
sphere_tree=cap_tree(sphere_tree);
cyl_tree.D=cyl_tree.D+1;
cyl_tree=cap_tree(cyl_tree);
cone_tree.D=cone_tree.D+1;
cone_tree=cap_tree(cone_tree);

Displacements=[0 , 0 , 0 ;  300 , 000 , 0 ; 600 , 0 , -150];

figure
hold on
hp= plot_tree(sphere_tree, [1 1 1], Displacements(1,:), [], [], '-b1');
set(hp,'facecolor',colDilution(2)*Colours(1,:)/256,'edgecolor','none');
hp= plot_tree(cyl_tree, [1 1 1], Displacements(2,:), [], [], '-b1');
set(hp,'facecolor',colDilution(3)*Colours(1,:)/256,'edgecolor','none');
hp= plot_tree(cone_tree, [1 1 1], Displacements(3,:), [], [], '-b1');
set(hp,'facecolor',colDilution(4)*Colours(1,:)/256,'edgecolor','none');

view([8.10000012278316 2.58619349745366])
axis equal off tight

saveName = '..\panels\fig_2\fig_2a_Shape_Morphs';

set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[9 4]);
tprint(saveName,'-SHR -eps',[9 4]);
tprint(saveName,'-SHR -tif',[9 4]);

%% 
%%%%%%%%%%%%%%%%%%%%% Panel A2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shape vec estimates
clear
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/Colours.mat')
load('../data/fig_2/shape_syns.mat')

colDilution=[1 , 0.8, 0.6 , 0.4];
res=50;
figure
hold on
plot(linspace(0,10),linspace(0,10),'black','LineWidth',0.5)
 
[PredGrid,MeanGrid,seGrid]=collate_synEst(shape_vec(1,:,:));
errorbar(PredGrid,MeanGrid,seGrid,'MarkerFaceColor','none','MarkerEdgeColor',colDilution(1)*Colours(1,:)/256,'Color',colDilution(1)*Colours(1,:)/256,'Marker','square','MarkerSize',1,'LineStyle','none','LineWidth',0.5,'LineStyle','none','CapSize',1)
[PredGrid,MeanGrid,seGrid]=collate_synEst(shape_vec(2,:,:));
errorbar(PredGrid,MeanGrid,seGrid,'MarkerFaceColor','none','MarkerEdgeColor',colDilution(2)*Colours(1,:)/256,'Color',colDilution(2)*Colours(1,:)/256,'Marker','o','MarkerSize',1,'LineStyle','none','LineWidth',0.5,'LineStyle','none','CapSize',1)
[PredGrid,MeanGrid,seGrid]=collate_synEst(shape_vec(3,:,:));
errorbar(PredGrid,MeanGrid,seGrid,'MarkerFaceColor','none','MarkerEdgeColor',colDilution(3)*Colours(1,:)/256,'Color',colDilution(3)*Colours(1,:)/256,'Marker','diamond','MarkerSize',1,'LineStyle','none','LineWidth',0.5,'LineStyle','none','CapSize',1)
[PredGrid,MeanGrid,seGrid]=collate_synEst(shape_vec(4,:,:));
errorbar(PredGrid,MeanGrid,seGrid,'MarkerFaceColor','none','MarkerEdgeColor',colDilution(4)*Colours(1,:)/256,'Color',colDilution(4)*Colours(1,:)/256,'Marker','^','MarkerSize',1,'LineStyle','none','LineWidth',0.5,'LineStyle','none','CapSize',1)

ylim([0 12])
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,5,10],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_2\fig_2a_Shape_syns');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
 
figure
hold on
 ShortNames={'Cube','Sphere','Cylinder','Cone'};
Locs=[(4:-1:1)' , zeros(4,1) ];

hold on
scatter(Locs(1,2),Locs(1,1),12,colDilution(1)*Colours(1,:)/256,'Marker','square');
scatter(Locs(2,2),Locs(2,1),8,colDilution(2)*Colours(1,:)/256,'Marker','o');
scatter(Locs(3,2),Locs(3,1),8,colDilution(3)*Colours(1,:)/256,'Marker','diamond');
scatter(Locs(4,2),Locs(4,1),8,colDilution(4)*Colours(1,:)/256,'Marker','^');

for ind=1:4
    text(Locs(ind,2)+0.35,Locs(ind,1),strcat(ShortNames{ind}),'FontSize',7,'FontName','helvetica','Color',colDilution(ind)*Colours(1,:)/256)
end
xlim([-1 3])
ylim([-1.5 6.5])
axis off
 saveName = strcat('..\panels\fig_2\fig_2a_Legend');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[2 4]);
 tprint(saveName,'-SHR -eps',[2 4]);
 tprint(saveName,'-SHR -tif',[2 4]);

 
 %% %% 
%%%%%%%%%%%%%%%%%%%%% Panel B1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BF  morphologies
clear
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/Colours.mat')
load('../data/fig_2/morphs_fig_2b.mat')

colDilution=[1 , 0.8, 0.6 , 0.4];
Displacements=[0 , 0 , 0 ;  300 , 000 , 0];

figure
hold on
hp= plot_tree(cap_tree(bf1_tree), [1 1 1], Displacements(1,:), [], [], '-b1');
set(hp,'facecolor',colDilution(1)*Colours(2,:)/256,'edgecolor','none');
hp= plot_tree(cap_tree(bf2_tree), [1 1 1], Displacements(2,:), [], [], '-b1');
set(hp,'facecolor',colDilution(4)*Colours(2,:)/256,'edgecolor','none');

view([8.10000012278316 2.58619349745366])
axis equal off tight

saveName = '..\panels\fig_2\fig_2b_BF_Morphs';

set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[6 4]);
tprint(saveName,'-SHR -eps',[6 4]);
tprint(saveName,'-SHR -tif',[6 4]);

bf1_tree=cap_tree(bf1_tree);
bf2_tree=cap_tree(bf2_tree);

%%  %% %% 
%%%%%%%%%%%%%%%%%%%%% Panel B2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BF Sholls and synaptic contacts
clear
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/Colours.mat')
load('../data/fig_2/bf_syns.mat')

res=50;
colDilution=[1 , 0.8, 0.6 , 0.4];
dd=linspace(0,400,res);
figure
hold on
for bf_ind=1:4
   plot(dd,ShollVecs(:,bf_ind),'Color',colDilution(bf_ind)*Colours(2,:)/256,'LineWidth',0.5);
end
set(gca,'ActivePositionProperty','position','XTick',[0,200,400],'YTick',[0,25,50],'XTickLabel',[0,100,200],'YTickLabel',[0,25,50],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_2\fig_2b_bf_sholls');   
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);
 

figure
hold on
plot(linspace(0,10),linspace(0,10),'black','LineWidth',0.5)
for bf_ind=1:4
    [PredGrid,MeanGrid,seGrid]=collate_synEst(bf_vec(bf_ind,:,:));
    errorbar(PredGrid,MeanGrid,seGrid,'o','MarkerFaceColor','none','MarkerSize',1,'MarkerEdgeColor',colDilution(bf_ind)*Colours(2,:)/256,'Color',colDilution(bf_ind)*Colours(2,:)/256,'Color',colDilution(bf_ind)*Colours(2,:)/256,'LineWidth',0.5,'LineStyle','none','CapSize',1);
end
ylim([0,12])
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,5,10],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_2\fig_2b_bf_syns');   
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

figure
hold on
ShortNames={'bf=0.2','bf=0.4','bf=0.6','bf=0.8'};
Locs=[(4:-1:1)' , zeros(4,1) ];

for bf_ind=1:4
    scatter(Locs(bf_ind,2),Locs(bf_ind,1),8,colDilution(bf_ind)*Colours(2,:)/256,'filled');
end

for ind=1:4
    text(Locs(ind,2)+0.35,Locs(ind,1),strcat(ShortNames{ind}),'FontSize',7,'FontName','helvetica','Color',colDilution(ind)*Colours(2,:)/256)
end
xlim([-1 3])
ylim([-1.5 6.5])

axis off
saveName = strcat('..\panels\fig_2\fig_2b_Legend');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[2 4]);
tprint(saveName,'-SHR -eps',[2 4]);
tprint(saveName,'-SHR -tif',[2 4]);

 %% 
%%%%%%%%%%%%%%%%%%%%% Panel C1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soma separation morphologies
clear
addpath('../functions/trees')
addpath('../functions')
load('../data/fig_2/morphs_fig_2c.mat')
load('../data/Colours.mat')

colDilution=[1 , 0.8, 0.6 , 0.4];
dend_tree.D=dend_tree.D+1;

figure
hold on
hp = plot_tree (cap_tree(dend_tree),[1 1 1], [0 0 000], [], [], '-b1');
set(hp,'facecolor',Colours(3,:)/256,'edgecolor','none');
hp= plot_tree (cap_tree(tran_tree(axonal_tree,[100 0 0])),[1 1 1], [0 0 000], [], [], '-b1');
set(hp,'facecolor',[0 0 0],'edgecolor', 'none');
axis equal off tight
saveName = '..\panels\fig_2\fig_2c_Intersoma_Morphs';
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[6 4]);
tprint(saveName,'-SHR -eps',[6 4]);
tprint(saveName,'-SHR -tif',[6 4]);


figure
hold on
xV=linspace(0,50,100);
yV=ones(100,1);
z = zeros(size(xV));
col=zeros(2,length(xV),3);
for pos_ind=1:length(xV)
   col(1,pos_ind,:)=(1.2-xV(pos_ind)/40)*Colours(3,:)/256;
   col(2,pos_ind,:)=(1.2-xV(pos_ind)/40)*Colours(3,:)/256; 
end
surface([xV;xV],[yV';yV'],[z;z],col,'facecol','no','edgecol','interp','linew',6);
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,5,10],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
axis off
saveName = strcat('..\panels\fig_2\fig_2c_Morphs_bar');   
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

%%  %% 
%%%%%%%%%%%%%%%%%%%%% Panel C2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soma separation synaptic contacts
clear
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/Colours.mat')
load('../data/fig_2/sep_syns.mat')

colDilution=[1 , 0.8, 0.6 , 0.4];
soma_seps=[ 10 , 20 , 30 , 40];
nsep=length(soma_seps);
Ests=zeros(nsep,2);
for sep_ind=1:nsep
   iVec=soma_vec(sep_ind,:,1);
   iVec(isnan(iVec))=0;
   Ests(sep_ind,1)=mean(iVec(:));
   Ests(sep_ind,2)=sqrt(var(iVec(:))/length(iVec));
end
figure
hold on
for sep_ind=1:4
    errorbar(soma_seps(sep_ind),Ests(sep_ind,1),Ests(sep_ind,2),'o','MarkerSize',3,'MarkerFaceColor',colDilution(sep_ind)*Colours(3,:)/256,'MarkerEdgeColor',colDilution(sep_ind)*Colours(3,:)/256,'Color',colDilution(sep_ind)*Colours(3,:)/256,'Color',colDilution(sep_ind)*Colours(3,:)/256,'LineWidth',0.5,'LineStyle','none','CapSize',1);
end
xV=linspace(0,50,100);
yV=interp1(soma_seps(:),2/pi*Ests(:,1),xV,'spline','extrap');
z = zeros(size(xV));

col=zeros(2,length(xV),3);
for pos_ind=1:length(xV)
    col(1,pos_ind,:)=(1.2-xV(pos_ind)/40)*Colours(3,:)/256;
   col(2,pos_ind,:)=(1.2-xV(pos_ind)/40)*Colours(3,:)/256; 
end
surface([xV;xV],[yV;yV],[z;z],col,'facecol','no','edgecol','interp','linew',2);


xlim([5 45])
ylim([3 5])
set(gca,'ActivePositionProperty','position','XTick',0:10:50,'YTick',[3,4,5],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_2\fig_2c_sep_Mean');   
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);


figure
hold on
plot(linspace(0,10),linspace(0,10),'black','LineWidth',0.5)
for sep_ind=1:4
    [PredGrid,MeanGrid,seGrid]=collate_synEst(soma_vec(sep_ind,:,:));
    errorbar(PredGrid,MeanGrid,seGrid,'o','MarkerFaceColor','none','MarkerSize',1,'MarkerEdgeColor',colDilution(sep_ind)*Colours(3,:)/256,'Color',colDilution(sep_ind)*Colours(3,:)/256,'Color',colDilution(sep_ind)*Colours(3,:)/256,'LineWidth',0.5,'LineStyle','none','CapSize',1);
end
ylim([0 12])
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,5,10],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_2\fig_2c_sep_syns');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

figure
hold on
ShortNames={'10\mum','20\mum','30\mum','40\mum'};
Locs=[(4:-1:1)' , zeros(4,1) ];
for bf_ind=1:4
    scatter(Locs(bf_ind,2),Locs(bf_ind,1),8,colDilution(bf_ind)*Colours(3,:)/256,'filled');
end
for ind=1:4
    text(Locs(ind,2)+0.35,Locs(ind,1),strcat(ShortNames{ind}),'FontSize',7,'FontName','helvetica','Color',colDilution(ind)*Colours(3,:)/256)
end
xlim([-1 3])
ylim([-1.5 6.5])
axis off
saveName = strcat('..\panels\fig_2\fig_2c_Legend');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[2 4]);
tprint(saveName,'-SHR -eps',[2 4]);
tprint(saveName,'-SHR -tif',[2 4]);

%% 
%%%%  %% 
%%%%%%%%%%%%%%%%%%%%% Panel D1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSCAM null morphs
clear
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/Colours.mat')
load('../data/fig_2/morphs_fig_2d.mat')

col1=Colours(4,:)/256;
colDilution=[1 , 0.8, 0.6 , 0.4];

for DS_ind=1:4
    DSTrees{DS_ind}=cap_tree(DSTrees{DS_ind});
end
figure
hp= plot_tree(DSTrees{1}, [1 1 1], [0 0 000], [], [], '-b1');
set(hp,'facecolor',colDilution(1)*col1,'edgecolor','none');
hp= plot_tree(tran_tree(DSTrees{2},[200 0 0]), [1 1 1], [0 0 000], [], [], '-b1');
set(hp,'facecolor',colDilution(2)*col1,'edgecolor','none');
hp= plot_tree(tran_tree(DSTrees{3},[0 -200 0]), [1 1 1], [0 0 000], [], [], '-b1');
set(hp,'facecolor',colDilution(3)*col1,'edgecolor','none');
hp= plot_tree(tran_tree(DSTrees{4},[200 -200 0]), [1 1 1], [0 0 000], [], [], '-b1');
set(hp,'facecolor',colDilution(4)*col1,'edgecolor','none');
view(2);
axis off
axis equal
saveName = '..\panels\fig_2\fig_2d_DSCAM_Morphs';
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[8 8]);
tprint(saveName,'-SHR -eps',[8 8]);
tprint(saveName,'-SHR -tif',[8 8]);

%%
%%%%  %% 
%%%%%%%%%%%%%%%%%%%%% Panel D2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSCAM null
clear
load('../data/Colours.mat')
load('../data/fig_2/DSCAM_syns.mat')


colDilution=[1 , 0.8, 0.6 , 0.4];
sz=1;
figure
hold on
nSamp=size(Len_vol_rat,1);

zVec=zeros(nSamp,4);
zVec(:,1)=4*(0:(nSamp-1))+1;
zVec(:,2)=4*(0:(nSamp-1))+2;
zVec(:,3)=4*(0:(nSamp-1))+3;
zVec(:,4)=4*(0:(nSamp-1))+4;

for DS_ind=1:4
    scatter3(Len_vol_rat(:,DS_ind,1),Len_vol_rat(:,DS_ind,2),zVec(:,DS_ind),sz,colDilution(DS_ind)*Colours(4,:)/256,'filled')%,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9);
end
xlim([0 , 5000]);
ylim([0 , 6*10^6]);
set(gca,'ActivePositionProperty','position','XTick',0:2500:5000,'YTick',[0,3,6]*10^6,'XTickLabel',0:2.5:5,'YTickLabel',[0,3,6],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');

saveName = strcat('..\panels\fig_2\fig_2d_LVrats');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
%plot(linspace(0,10),linspace(0,10),'black','LineWidth',0.5)
for DS_ind=1:4
    [PredGrid,MeanGrid,seGrid,~,~,~,nEach]=collate_synEst(DSCAM_vec(:,DS_ind,[1 , 4]),max(DSCAM_vec(:,DS_ind,1)));
    errorbar(PredGrid(nEach>10),MeanGrid(nEach>10),seGrid(nEach>10),'o','MarkerFaceColor','none','MarkerSize',1,'MarkerEdgeColor',colDilution(DS_ind)*Colours(4,:)/256,'Color',colDilution(DS_ind)*Colours(4,:)/256,'Color',colDilution(DS_ind)*Colours(4,:)/256,'LineWidth',0.5,'LineStyle','--','CapSize',1);
end
ylim([0,12])
xlim([0 , 5000]);
set(gca,'ActivePositionProperty','position','XTick',[0,2500,5000],'YTick',[0,5,10],'XTickLabel',[0,2.5,5],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_2\fig_2d_DSCAM_lengthCon');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
%plot(linspace(0,10),linspace(0,10),'black','LineWidth',0.5)
for DS_ind=1:4
    [PredGrid,MeanGrid,seGrid,~,~,~,nEach]=collate_synEst(DSCAM_vec(:,DS_ind,[2 , 4]),max(DSCAM_vec(:,DS_ind,2)));
    errorbar(PredGrid(nEach>10),MeanGrid(nEach>10),seGrid(nEach>10),'o','MarkerFaceColor','none','MarkerSize',1,'MarkerEdgeColor',colDilution(DS_ind)*Colours(4,:)/256,'Color',colDilution(DS_ind)*Colours(4,:)/256,'Color',colDilution(DS_ind)*Colours(4,:)/256,'LineWidth',0.5,'LineStyle','--','CapSize',1);
end
ylim([0,12])
xlim([1000000,6000000])
set(gca,'ActivePositionProperty','position','XTick',[1000000,3500000,6000000],'YTick',[0,5,10],'XTickLabel',[1,3.5,6],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_2\fig_2d_DSCAM_volCon');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
plot(linspace(0,10),linspace(0,10),'black','LineWidth',0.5)
for DS_ind=1:4
    [PredGrid,MeanGrid,seGrid]=collate_synEst(DSCAM_vec(:,DS_ind,[3 4]));
    errorbar(PredGrid,MeanGrid,seGrid,'o','MarkerFaceColor','none','MarkerSize',1,'MarkerEdgeColor',colDilution(DS_ind)*Colours(4,:)/256,'Color',colDilution(DS_ind)*Colours(4,:)/256,'Color',colDilution(DS_ind)*Colours(4,:)/256,'LineWidth',0.5,'LineStyle','none','CapSize',1);
end
ylim([0,12])
set(gca,'ActivePositionProperty','position','XTick',[0,5,10],'YTick',[0,5,10],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
saveName = strcat('..\panels\fig_2\fig_2d_DSCAM_syns');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

figure
hold on
ShortNames={'0','100','400','700'};
Locs=[(4:-1:1)' , zeros(4,1) ];

for DS_ind=1:4
    scatter(Locs( DS_ind,2),Locs( DS_ind,1),8,colDilution( DS_ind)*Colours(4,:)/256,'filled');
end

for ind=1:4
    text(Locs(ind,2)+0.35,Locs(ind,1),strcat(ShortNames{ind}),'FontSize',7,'FontName','helvetica','Color',colDilution(ind)*Colours(4,:)/256)
end
xlim([-1 3])
ylim([-1.5 6.5])

axis off
saveName = strcat('..\panels\fig_2\fig_2d_Legend');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[2 4]);
tprint(saveName,'-SHR -eps',[2 4]);
tprint(saveName,'-SHR -tif',[2 4]);


