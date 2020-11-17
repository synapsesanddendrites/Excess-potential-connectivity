 %% Mouse morphologies
 clear
load('../data/Colours.mat')
load('../data/fig_5/fig_5amouseMorphs.mat')

col1=Colours(1,:)/256;
col2=Colours(2,:)/256;

%[axcell,dencell]=cell_pair_plot(oMouse);

tree1=dencell;
tree2=axcell;

d1=tree1.D;
d2=tree2.D;

d1M=mean(d1(:));
d2M=mean(d2(:));

d1(d1>(3*d1M))=d1M;
d2(d2>(3*d1M))=d2M;

tree1.D=d1+2;
tree2.D=d2+1;

tree1.Y=-tree1.Y;
tree2.Y=-tree2.Y;

tree1=soma_tree(tree1,20);
tree2=soma_tree(tree2,20);

figure
hp= plot_tree(tree1, [1 1 1], [0 0 000], [], [], '-b1');
set(hp,'facecolor',col1,'edgecolor','none');
view(2);

hp= plot_tree(tree2, [1 1 1], [0 0 000], [], [], '-b1');
set(hp,'facecolor',col2,'edgecolor','none');
view(2);
hp= line ([50 100], [-100 -100]);
set(hp,'color',[0 0 0],'linewidth',1);

axis equal off tight

saveName = '..\panels\fig_5\fig_5a_Mouse_Example';

set(gcf,'renderer','painter','PaperPositionMode','manual');

tprint(saveName,'-SHR -jpg',[7 12]);
set(gcf,'renderer','painter');
tprint(saveName,'-SHR -eps',[7 12]);
tprint(saveName,'-SHR -tif',[7 12]);


%% Mouse grids
 clear
load('../data/Colours.mat')
load('../data/fig_5/fig_5mouse.mat')

maxCols=[4000000,16000,4000,120];
[parse_mouse,len_mouse]=parse_gridm(mouseGrid,maxCols);
xLabels={0,[],[],[],200,[],[],[],400,[],[],[],600,[],[],[],800};
yLabels{1}={0,[],[],[],1,[],[],[],2,[],[],[],3,[],[],[],4};
yLabels{2}={0,[],[],[],4,[],[],[],8,[],[],[],12,[],[],[],16};
yLabels{3}={0,[],[],[],1,[],[],[],2,[],[],[],3,[],[],[],4};
yLabels{4}={0,[],[],[],30,[],[],[],60,[],[],[],90,[],[],[],120};
names={'V','La','Ld','N'};
for col_ind=1:4
    figure
    hold on
    thisMat=parse_mouse{col_ind,1};
    allVals=thisMat(:);
    allVals=sort(allVals);
    allVals=unique(allVals);
    thisGrid=parse_mouse{col_ind,2};
    for ind_sep=1:length(len_mouse)
        for ind_Meas=1:length(len_mouse)
            if thisMat(ind_sep,ind_Meas)>0
                Gsize=interp1(allVals,linspace(0,1,length(allVals)),thisMat(ind_sep,ind_Meas));
                scatter3(len_mouse(ind_sep),thisGrid(ind_Meas),0,32*Gsize,(thisMat(ind_sep,ind_Meas)),'filled','square');
            end
        end
    end
    colMean=parse_mouse{col_ind,3};
    smoothBase=linspace(0,800,100);
    smoothMeans=interp1(len_mouse,colMean,smoothBase,'spline');
    
    plot3(smoothBase,smoothMeans,ones(100,1),'Color','black','LineWidth',0.5);
    
    
    colormap('Jet')
    xlim([0 800])
    ylim([0 maxCols(col_ind)])
    caxis([0 1])
    set(gca,'ActivePositionProperty','position','XTick',linspace(0,800,length(len_mouse)+1),'YTick',linspace(0,maxCols(col_ind),length(len_mouse)+1),'XTickLabels',xLabels,'YTickLabels',yLabels{col_ind},'ticklength',[0.02 0.02],'fontsize',8,'fontname','helvetica','XGrid','on','YGrid','on','tickdir','out');
    
    saveName = strcat('..\panels\fig_5\fig_5b_',num2str(col_ind),names{col_ind});
    set(gcf,'renderer','painter','PaperPositionMode','manual');
    tprint(saveName,'-SHR -jpg',[4 4]);
    tprint(saveName,'-SHR -eps',[4 4]);
    tprint(saveName,'-SHR -tif',[4 4]);       
end

cbh=colorbar('h');
caxis([0 1])
set(cbh,'YTick',0:0.5:1,'fontsize',8,'fontname','helvetica')

saveName = strcat('..\panels\fig_5\fig_5b_ColorBar');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

    
    %% %% Human morphologies
    
clear
load('../data/Colours.mat')
load('../data/fig_5/fig_5ahumanMorphs.mat')

col1=Colours(1,:)/256;
col2=Colours(2,:)/256;

tree1=dencell;
tree2=axcell;

d1=tree1.D;
d2=tree2.D;

d1M=mean(d1(:));
d2M=mean(d2(:));

d1(d1>(3*d1M))=d1M;
d2(d2>(3*d1M))=d2M;

tree1.D=d1+1;
tree2.D=d2+2;

tree1.Y=-tree1.Y;
tree2.Y=-tree2.Y;

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

axis equal off tight

saveName = '..\panels\fig_5\fig_5a_Human_Example';
set(gcf,'renderer','painter','PaperPositionMode','manual');

tprint(saveName,'-SHR -jpg',[7 12]);
set(gcf,'renderer','painter');
tprint(saveName,'-SHR -eps',[7 12]);
tprint(saveName,'-SHR -tif',[7 12]);
 %% Human heat maps
 clear
load('../data/Colours.mat')
load('../data/fig_5/fig_5human.mat')

maxCols=[30000000,25000,16000,80];
[parse_human,len_human]=parse_gridh(HumanGrid,maxCols);
xLabels={0,[],[],[],400,[],[],[],800,[],[],[],1200,[],[],[],1600};
yLabels{1}={0,[],[],[],7.5,[],[],[],15,[],[],[],22.5,[],[],[],30};
yLabels{2}={0,[],[],[],6.25,[],[],[],12.5,[],[],[],18.75,[],[],[],25};
yLabels{3}={0,[],[],[],4,[],[],[],8,[],[],[],12,[],[],[],16};
yLabels{4}={0,[],[],[],20,[],[],[],40,[],[],[],60,[],[],[],80};
names={'V','La','Ld','N'};
for col_ind=1:4
    figure
    hold on
    thisMat=parse_human{col_ind,1};
    allVals=thisMat(:);
    allVals=sort(allVals);
    allVals=unique(allVals);
    thisGrid=parse_human{col_ind,2};
    for ind_sep=1:length(len_human)
        for ind_Meas=1:length(len_human)
            if thisMat(ind_sep,ind_Meas)>0
                Gsize=interp1(allVals,linspace(0,1,length(allVals)),thisMat(ind_sep,ind_Meas));
                scatter3(len_human(ind_sep),thisGrid(ind_Meas),0,32*Gsize,(thisMat(ind_sep,ind_Meas)),'filled','square');
            end
        end
    end
    colMean=parse_human{col_ind,3};
    smoothBase=linspace(0,1600,100);
    smoothMeans=interp1(len_human,colMean,smoothBase,'spline');
    
    plot3(smoothBase,smoothMeans,ones(100,1),'Color','black','LineWidth',0.5);
    
    colormap('Jet')
    xlim([0 1200])
    ylim([0 maxCols(col_ind)])
    caxis([0 1])
    set(gca,'ActivePositionProperty','position','XTick',linspace(0,1200,length(len_human)+1),'YTick',linspace(0,maxCols(col_ind),length(len_human)+1),'XTickLabels',xLabels,'YTickLabels',yLabels{col_ind},'ticklength',[0.02 0.02],'fontsize',8,'fontname','helvetica','XGrid','on','YGrid','on','tickdir','out');
    
    saveName = strcat('..\panels\fig_5\fig_5d_',num2str(col_ind),names{col_ind});
    set(gcf,'renderer','painter','PaperPositionMode','manual');
    tprint(saveName,'-SHR -jpg',[4 4]);
    tprint(saveName,'-SHR -eps',[4 4]);
    tprint(saveName,'-SHR -tif',[4 4]);    
end

    