load('../data/fig_5/supdata.mat') 
load('../data/Colours.mat')


sz=6;

figure
hold on


scatter(MouseAx(:,1),MouseAx(:,2),sz,Colours(2,:)/256,'filled');
xlim([0 , 40000]);
ylim([0 , 3*10^7]);
set(gca,'ActivePositionProperty','position','XTick',0:20000:40000,'YTick',[0,1.5,3]*10^7,'XTickLabel',0:20:40,'YTickLabel',[0,1.5,3],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');


saveName = strcat('..\panels\fig_5\fig_S5_MouseAxonsLV');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
  tprint(saveName,'-SHR -tif',[4 4]);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on

scatter(MouseDen(:,1),MouseDen(:,2),sz,Colours(1,:)/256,'filled');
xlim([0 , 10000]);
ylim([0 , 10^7]);
set(gca,'ActivePositionProperty','position','XTick',0:5000:10000,'YTick',[0,0.5,1]*10^7,'XTickLabel',0:5:10,'YTickLabel',[0,0.5,1],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');

saveName = strcat('..\panels\fig_5\fig_S5_MouseDensLV');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
  tprint(saveName,'-SHR -tif',[4 4]);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on

scatter(HumanAx(:,1),HumanAx(:,2),sz,Colours(2,:)/256,'filled');
xlim([0 , 50000]);
ylim([0 , 15*10^7]);
set(gca,'ActivePositionProperty','position','XTick',0:25000:50000,'YTick',[0,7.5,15]*10^7,'XTickLabel',0:25:50,'YTickLabel',[0,7.5,15],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');

saveName = strcat('..\panels\fig_5\fig_S5_HumanAxonsLV');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
  tprint(saveName,'-SHR -tif',[4 4]);
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on

scatter(HumanDen(:,1),HumanDen(:,2),sz,Colours(1,:)/256,'filled');
 xlim([0 , 20000]);
 ylim([0 , 5*10^7]);
set(gca,'ActivePositionProperty','position','XTick',0:10000:20000,'YTick',[0,2.5,5]*10^7,'XTickLabel',0:10:20,'YTickLabel',[0,2.5,5],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');


saveName = strcat('..\panels\fig_5\fig_S5_HumanDensLV');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
  tprint(saveName,'-SHR -tif',[4 4]);
  
 %% 
clear
load('../data/fig_5/supdata.mat') 
load('../data/Colours.mat')
res=100;

nMouse=length(MouseDens);
mouseAxall=[];
mouseDenall=[];
for mouse_ind=1:nMouse
    mouseAxall=[mouseAxall ; MouseDens{mouse_ind,2}];
    mouseDenall=[mouseDenall ; MouseDens{mouse_ind,1}];
end
maxALL=max(max(mouseAxall),max(mouseDenall));
mouseGrid=linspace(0,maxALL,res+1);
slicesize=mouseGrid(2);

[NmouseAx,MouseAXedges] = histcounts(mouseAxall,mouseGrid);
[NmouseDen,MouseDENedges] = histcounts(mouseDenall,mouseGrid);
Xmouse=(mouseGrid(1:res)+mouseGrid(2:end))/2;

figure
hold on

plot(Xmouse,NmouseAx/slicesize,'Color',Colours(2,:)/256,'LineWidth',2)
plot(Xmouse,NmouseDen/slicesize,'Color',Colours(1,:)/256,'LineWidth',2)
xlim([0 , 1200]);
ylim([0 , 4000]);
set(gca,'ActivePositionProperty','position','XTick',0:400:1200,'YTick',0:2000:4000,'XTickLabel',0:0.4:1.2,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');

saveName = strcat('..\panels\fig_5\fig_S5_MouseNeurite');   
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[NmouseSoma] = histcounts(MouseSoma,mouseGrid);

figure
hold on

stem(Xmouse,NmouseSoma/slicesize,'filled','Color',[0 0 0],'LineWidth',0.5,'MarkerSize',1)

xlim([0 , 1200]);
ylim([0 , 0.75]);
set(gca,'ActivePositionProperty','position','XTick',0:400:1200,'YTick',0:0.25:0.75,'XTickLabel',0:0.4:1.2,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');

saveName = strcat('..\panels\fig_5\fig_S5_MouseSoma');   
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

 %% 
clear
load('../data/fig_5/supdata.mat') 
load('../data/Colours.mat')
res=100;

nHuman=length(HumanDens);
HumanAxall=[];
HumanDenall=[];
for Human_ind=1:nHuman
    HumanAxall=[HumanAxall ; HumanDens{Human_ind,2}];
    HumanDenall=[HumanDenall ; HumanDens{Human_ind,1}];
end
maxALL=max(max(HumanAxall),max(HumanDenall));
HumanGrid=linspace(0,maxALL,res+1);
slicesize=HumanGrid(2);

[NHumanAx,HumanAXedges] = histcounts(HumanAxall,HumanGrid);
[NHumanDen,HumanDENedges] = histcounts(HumanDenall,HumanGrid);
XHuman=(HumanGrid(1:res)+HumanGrid(2:end))/2;

figure
hold on

plot(XHuman,NHumanAx/slicesize,'Color',Colours(2,:)/256,'LineWidth',2)
plot(XHuman,NHumanDen/slicesize,'Color',Colours(1,:)/256,'LineWidth',2)
xlim([0 , 2500]);
ylim([0 , 2000]);
set(gca,'ActivePositionProperty','position','XTick',0:1250:2500,'YTick',0:1000:2000,'XTickLabel',0:1.25:2.5,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');

saveName = strcat('..\panels\fig_5\fig_S5_HumanNeurite');   
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[NHumanSoma] = histcounts(HumanSoma,HumanGrid);

figure
hold on

stem(XHuman,NHumanSoma/slicesize,'filled','Color',[0 0 0],'LineWidth',0.5,'MarkerSize',1)

xlim([0 , 2500]);
ylim([0 , 0.3]);
set(gca,'ActivePositionProperty','position','XTick',0:1250:2500,'YTick',0:0.1:0.3,'XTickLabel',0:1.25:2.5,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');

saveName = strcat('..\panels\fig_5\fig_S5_HumanSoma');   
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);