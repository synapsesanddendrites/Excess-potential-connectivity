%%%%%%%%%%%%%%%%%%%%% Panel B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shared volume schematic
clear
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/Colours.mat')
load('../data/fig_1/morphs_fig_1.mat')

col1=Colours(1,:)/256;
col2=Colours(2,:)/256;

bndcol=[0.75 0.75 0.75];
[c1]=convexity_tree(dend_tree);
[c2]=convexity_tree(axonal_tree);

bound1 = boundary_tree(dend_tree,'-3d',c1); % Get boundary
bound2 = boundary_tree(axonal_tree,'-3d',c2); % Get boundary

tree1=resample_tree(dend_tree,1);
tree2=resample_tree(axonal_tree,1);

x1=tree1.X;
y1=tree1.Y;
z1=tree1.Z;
x2=tree2.X;
y2=tree2.Y;
z2=tree2.Z;

Is1in2=intriangulation(bound2.Vertices,bound2.Faces,[x1,y1,z1]);
Is2in1=intriangulation(bound1.Vertices,bound1.Faces,[x2,y2,z2]);

Len1=nnz(Is1in2)*mean(len_tree(tree1));
Len2=nnz(Is2in1)*mean(len_tree(tree2));

X=[x1(Is1in2); x2(Is2in1)];
Y=[y1(Is1in2); y2(Is2in1)];
Z=[z1(Is1in2); z2(Is2in1)];

[k,sharedVol]=boundary(X,Y,Z,1-(c1+c2)/2);

figure
F=gcf;
h=trisurf(k,X,Y,Z);
[rh]=reducepatch(h,0.5);
rh.Vertices=rh.vertices;
rh.Faces=rh.faces;
pbound=rh;
close(F)
%
allowedAx=intriangulation(pbound.vertices,pbound.faces,[axonal_tree.X,axonal_tree.Y,axonal_tree.Z]);
allowedDen=intriangulation(pbound.vertices,pbound.faces,[dend_tree.X,dend_tree.Y,dend_tree.Z]);

[ibound,~]=boundary(X,Y,Z,1-(c1+c2)/2);
figure
hold on
fv = struct('vertices',[X,Y,Z],'faces',ibound);
h2= trisurf(ibound,X,Y,Z,'FaceColor',bndcol,'EdgeColor','none','FaceAlpha',0.2);


tree1.D=tree1.D+1;
tree2.D=tree2.D+1;


hp = plot_tree (tree1,[1 1 1], [0 0 000], Is1in2, 16, '-b1');
set(hp,'facecolor',col1,'edgecolor','none');
hp= plot_tree (tree2,[1 1 1], [0 0 000], Is2in1, 16, '-b1');
set(hp,'facecolor',col2,'edgecolor', 'none');
hp = plot_tree (tree1,[1 1 1], [0 0 000], (1-Is1in2)==1, 16, '-b1');
set(hp,'facecolor',col1,'edgecolor','none','FaceAlpha',0.2);
hp= plot_tree (tree2,[1 1 1], [0 0 000], (1-Is2in1)==1, 16, '-b1');
set(hp,'facecolor',col2,'edgecolor', 'none','FaceAlpha',0.2);

view(2);
axis equal off tight
set(gcf,'renderer', 'opengl');

savename='../panels/fig_1/fig_1b_Shared_volume';
tprint(savename,'-SHR -jpg',[6 6]);
set(gcf,'renderer','painter');
tprint(savename,'-SHR -eps',[6 6]);
tprint(savename,'-SHR -tif',[6 6]);

%% 
%%%%%%%%%%%%%%%%%%%%% Panel C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Trend in terms of dendrite length
clear
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/Colours.mat')
load('../data/fig_1/dendrite_panel')

lD=2.4e+03;
lA=3.0e+03;
vS=2.4e+06;
spineVec=[1.5,2.5,3.5,4.5];


nStep=10;

%%%%%%%%%%%%% Sort values %%%%%%%%%%%%%%%%%%%%%%
minDenLen=min(denVals(:,2));
maxDenLen=4000;

DenGrid=linspace(minDenLen,maxDenLen,nStep+1);
DenMids=(DenGrid(2:(nStep+1))+DenGrid(1:nStep))/2;

theseINDsCell=cell(nStep,1);
for stpInd=1:nStep
    DenMin=DenGrid(stpInd);
    DenMax=DenGrid(stpInd+1);
    
    theseINDS=denVals(:,2)>=DenMin & denVals(:,2)<=DenMax;
    nVals=nnz(theseINDS);
    theseINDsCell{stpInd}=find(theseINDS);    
end

allNums=cell(nStep,length(spineVec));
allMeans=zeros(nStep,length(spineVec));
allSEs=zeros(nStep,length(spineVec));
allSDs=zeros(nStep,length(spineVec));
allSDSEs=zeros(nStep,length(spineVec));
for spInd=1:length(spineVec)
    for stpInd=1:nStep  
        theseINDS=theseINDsCell{stpInd};
        theseSyns=denTrend(theseINDS,2,spInd);
        
        theseSyns=theseSyns(:);
        allNums{stpInd,spInd}=theseSyns;
        
        allMeans(stpInd,spInd)=mean(theseSyns);
        allSEs(stpInd,spInd)=sqrt(var(theseSyns)/length(theseINDS));
        
        thisSD=sqrt(var(theseSyns));
        allSDs(stpInd,spInd)=thisSD;
        allSDSEs(stpInd,spInd)=thisSD*sqrt(2/(length(theseINDS)-1));
    end   
end

%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%
figure
hold on
iX=linspace(minDenLen,maxDenLen,100);
for spInd=1:length(spineVec)
    errorbar(DenMids,allMeans(:,spInd),allSEs(:,spInd),'MarkerFaceColor',Colours(spInd,:)/256,'MarkerEdgeColor',Colours(spInd,:)/256,'Marker','square','MarkerSize',1,'CapSize',1,'LineStyle','none');
    iY=pi/2*lA*iX*spineVec(spInd)/vS;
    plot(iX,iY,'Color',Colours(spInd,:)/256,'LineStyle','-','LineWidth',0.5)
end
xlim([0 4500])
ylim([0 50])
set(gca,'ActivePositionProperty','position','XTick',0:2000:4000,'XTickLabel',0:2:4,'YTick',0:25:50,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
savename = '../panels/fig_1/fig_1c_Dendrite_trend';
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(savename,'-SHR -jpg',[4 4]);
tprint(savename,'-SHR -eps',[4 4]);
tprint(savename,'-SHR -tif',[4 4]);

figure
hold on
for spInd=1:length(spineVec)
    errorbar(DenMids,allSDs(:,spInd),allSDSEs(:,spInd),'MarkerFaceColor',Colours(spInd,:)/256,'MarkerEdgeColor',Colours(spInd,:)/256,'Marker','square','MarkerSize',1,'CapSize',1,'LineStyle','-','Color',Colours(spInd,:)/256);
end
xlim([0 4500])
ylim([0 50])
set(gca,'ActivePositionProperty','position','XTick',0:2000:4000,'XTickLabel',0:2:4,'YTick',0:25:50,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
savename = '../panels/fig_1/fig_1c_Dendrite_variance_trend';
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(savename,'-SHR -jpg',[4 4]);
tprint(savename,'-SHR -eps',[4 4]);
tprint(savename,'-SHR -tif',[4 4]);



%% 
%%%%%%%%%%%%%%%%%%%%% Panel D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trend in terms of axon length

clear
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/Colours.mat')
load('../data/fig_1/axon_panel.mat')

lD=2.4e+03;
lA=3.0e+03;
vS=2.4e+06;
spineVec=[1.5,2.5,3.5,4.5];

nStep=10;

%%%%%%%%%%%%% Sort values %%%%%%%%%%%%%%%%%%%%%%
minAxLen=1000;
maxAxLen=5000;

axGrid=linspace(minAxLen,maxAxLen,nStep+1);
axMids=(axGrid(2:(nStep+1))+axGrid(1:nStep))/2;

theseINDsCell=cell(nStep,1);
for stpInd=1:nStep
    axMin=axGrid(stpInd);
    axMax=axGrid(stpInd+1);
    
    theseINDS=axVals(:,3)>=axMin & axVals(:,3)<=axMax;
    nVals=nnz(theseINDS);
    theseINDsCell{stpInd}=find(theseINDS);    
end

allNums=cell(nStep,length(spineVec));
allMeans=zeros(nStep,length(spineVec));
allSEs=zeros(nStep,length(spineVec));
allSDs=zeros(nStep,length(spineVec));
allSDSEs=zeros(nStep,length(spineVec));

for spInd=1:length(spineVec)
    for stpInd=1:nStep  
        theseINDS=theseINDsCell{stpInd};
        theseSyns=axTrend(theseINDS,2,spInd);
        
        theseSyns=theseSyns(:);
        allNums{stpInd,spInd}=theseSyns;
        
        allMeans(stpInd,spInd)=mean(theseSyns);
        allSEs(stpInd,spInd)=sqrt(var(theseSyns)/length(theseINDS));
        
        thisSD=sqrt(var(theseSyns));
        allSDs(stpInd,spInd)=thisSD;
        allSDSEs(stpInd,spInd)=thisSD*sqrt(2/(length(theseINDS)-1));
    end   
end

%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%
figure
hold on
iX=linspace(minAxLen,maxAxLen,100);
for spInd=1:length(spineVec)
    errorbar(axMids,allMeans(:,spInd),allSEs(:,spInd),'MarkerFaceColor',Colours(spInd,:)/256,'MarkerEdgeColor',Colours(spInd,:)/256,'Marker','square','MarkerSize',1,'CapSize',1,'LineStyle','none');
    iY=pi/2*lD*iX*spineVec(spInd)/vS;
    plot(iX,iY,'Color',Colours(spInd,:)/256,'LineStyle','-','LineWidth',0.5) 
end
xlim([0000 5500])
ylim([0 50])
set(gca,'ActivePositionProperty','position','XTick',0:2500:5000,'XTickLabel',0:2.5:5,'YTick',0:25:50,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');

savename = strcat('../panels/fig_1/fig_1d_Axonal_trend');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(savename,'-SHR -jpg',[4 4]);
tprint(savename,'-SHR -eps',[4 4]);
tprint(savename,'-SHR -tif',[4 4]);

figure
hold on
for spInd=1:length(spineVec)
    errorbar(axMids,allSDs(:,spInd),allSDSEs(:,spInd),'MarkerFaceColor',Colours(spInd,:)/256,'MarkerEdgeColor',Colours(spInd,:)/256,'Marker','square','MarkerSize',1,'CapSize',1,'LineStyle','-','Color',Colours(spInd,:)/256);
end
xlim([0 5500])
ylim([0 50])
set(gca,'ActivePositionProperty','position','XTick',0:2500:5000,'XTickLabel',0:2.5:5,'YTick',0:25:50,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
savename = '../panels/fig_1/fig_1d_Axonal_variance_trend';
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(savename,'-SHR -jpg',[4 4]);
tprint(savename,'-SHR -eps',[4 4]);
tprint(savename,'-SHR -tif',[4 4]);



%% 
%
%%%%%%%%%%%%%%%%%%%%% Panel E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trend in terms of volume

clear
addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/Colours.mat')
load('../data/fig_1/volume_panel.mat')


lD=2.4e+03;
lA=3.0e+03;
vS=2.4e+06;
spineVec=[1.5,2.5,3.5,4.5];

nStep=10;

%%%%%%%%%%%%% Sort values %%%%%%%%%%%%%%%%%%%%%%
minVol=1900000; 
maxVol=2900000;

VolGrid=linspace(minVol,maxVol,nStep+1);
VolMids=(VolGrid(2:(nStep+1))+VolGrid(1:nStep))/2;

theseINDsCell=cell(nStep,1);
for stpInd=1:nStep
    VolMin=VolGrid(stpInd);
    VolMax=VolGrid(stpInd+1);
    
    theseINDS=shVals(:,1)>=VolMin & shVals(:,1)<=VolMax;
    nVals=nnz(theseINDS);
    theseINDsCell{stpInd}=find(theseINDS);    
end

allNums=cell(nStep,length(spineVec));
allMeans=zeros(nStep,length(spineVec));
allSEs=zeros(nStep,length(spineVec));
allSDs=zeros(nStep,length(spineVec));
allSDSEs=zeros(nStep,length(spineVec));

for spInd=1:length(spineVec)
    for stpInd=1:nStep  
        theseINDS=theseINDsCell{stpInd};
        theseSyns=shTrend(theseINDS,2,spInd);
        
        theseSyns=theseSyns(:);
        allNums{stpInd,spInd}=theseSyns;
        
        allMeans(stpInd,spInd)=mean(theseSyns);
        allSEs(stpInd,spInd)=sqrt(var(theseSyns)/length(theseINDS));
        
          
        thisSD=sqrt(var(theseSyns));
        allSDs(stpInd,spInd)=thisSD;
        allSDSEs(stpInd,spInd)=thisSD*sqrt(2/(length(theseINDS)-1));
    end   
end

%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%
figure
hold on
iX=linspace(minVol,maxVol,100);
for spInd=1:length(spineVec)
    errorbar(VolMids,allMeans(:,spInd),allSEs(:,spInd),'MarkerFaceColor',Colours(spInd,:)/256,'MarkerEdgeColor',Colours(spInd,:)/256,'Marker','square','MarkerSize',1,'CapSize',1,'LineStyle','none');
    iY=pi/2*lA*lD*spineVec(spInd)./iX;
    plot(iX,iY,'Color',Colours(spInd,:)/256,'LineStyle','-','LineWidth',0.5) 
end
xlim([1800000 3000000])
ylim([0 50])
set(gca,'ActivePositionProperty','position','XTick',1800000:600000:3000000,'XTickLabel',1.8:0.6:3,'YTick',0:25:50,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
savename = '../panels/fig_1/fig_1e_Volume_trend';
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(savename,'-SHR -jpg',[4 4]);
tprint(savename,'-SHR -eps',[4 4]);
tprint(savename,'-SHR -tif',[4 4]);

figure
hold on
for spInd=1:length(spineVec)
    errorbar(VolMids,allSDs(:,spInd),allSDSEs(:,spInd),'MarkerFaceColor',Colours(spInd,:)/256,'MarkerEdgeColor',Colours(spInd,:)/256,'Marker','square','MarkerSize',1,'CapSize',1,'LineStyle','-','Color',Colours(spInd,:)/256);
end
xlim([1800000 3000000])
ylim([0 50])
set(gca,'ActivePositionProperty','position','XTick',1800000:600000:3000000,'XTickLabel',1.8:0.6:3,'YTick',0:25:50,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');

savename = '../panels/fig_1/fig_1e_Volume_variance_trend';
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(savename,'-SHR -jpg',[4 4]);
tprint(savename,'-SHR -eps',[4 4]);
tprint(savename,'-SHR -tif',[4 4]);

 
 %% 
clear
load('../data/Colours.mat')

 ShortNames={'1.5','2.5','3.5','4.5'};
Locs=[(4:-1:1)' , zeros(4,1)];
figure
hold on
for ind=1:4
    scatter(Locs(ind,2),Locs(ind,1),12,Colours(ind,:)/256,'filled','Marker','square');
    text(Locs(ind,2)+0.35,Locs(ind,1),strcat(ShortNames{ind}),'FontSize',7,'FontName','helvetica')
end
xlim([-1 3])
ylim([0 6])
axis off
 savename ='../panels/fig_1/fig_1_Legend';   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(savename,'-SHR -jpg',[4 4]);
 tprint(savename,'-SHR -eps',[4 4]);
 tprint(savename,'-SHR -tif',[4 4]);