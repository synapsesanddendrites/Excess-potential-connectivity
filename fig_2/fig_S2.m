load('../data/fig_2/Shape_LV.mat')
load('../data/Colours.mat')


colDilution=[1 , 0.8, 0.6 , 0.4];
sz=6;
figure
hold on
nSamp=size(Len_vol_rat,2);
zVec1=4*(0:(nSamp-1))+1;
zVec2=4*(0:(nSamp-1))+2;
zVec3=4*(0:(nSamp-1))+3;
zVec4=4*(0:(nSamp-1))+4;

    scatter3(Len_vol_rat(1,:,1),Len_vol_rat(1,:,2),zVec1,sz,colDilution(1)*Colours(1,:)/256,'filled','Marker','square');
    scatter3(Len_vol_rat(2,:,1),Len_vol_rat(2,:,2),zVec2,sz,colDilution(2)*Colours(1,:)/256,'filled','Marker','o');
    scatter3(Len_vol_rat(3,:,1),Len_vol_rat(3,:,2),zVec3,sz,colDilution(3)*Colours(1,:)/256,'filled','Marker','diamond');
    scatter3(Len_vol_rat(4,:,1),Len_vol_rat(4,:,2),zVec4,sz,colDilution(4)*Colours(1,:)/256,'filled','Marker','^');


xlim([0 , 5000]);
ylim([0 , 10^7]);
set(gca,'ActivePositionProperty','position','XTick',0:2500:5000,'YTick',[0,5,10]*10^6,'XTickLabel',0:2.5:5,'YTickLabel',[0,5,10],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
%legend('Dendrite','Axon','boxoff','FontSize',6,'FontName','helvetica','Location','northwest');



saveName = strcat('..\panels\fig_2\fig_S2_aShapes');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
 
 
figure
hold on

 ShortNames={'Cube','Sphere','Cylinder','Cone'};
Locs=[(2:-1:1)' , zeros(2,1) ; (2:-1:1)' , 3*ones(2,1)];

hold on
scatter(Locs(1,2),Locs(1,1),8,colDilution(1)*Colours(1,:)/256,'filled','Marker','square');
scatter(Locs(2,2),Locs(2,1),8,colDilution(2)*Colours(1,:)/256,'filled','Marker','o');
scatter(Locs(3,2),Locs(3,1),8,colDilution(3)*Colours(1,:)/256,'filled','Marker','diamond');
scatter(Locs(4,2),Locs(4,1),8,colDilution(4)*Colours(1,:)/256,'filled','Marker','^');

for ind=1:4
    text(Locs(ind,2)+0.35,Locs(ind,1),strcat(ShortNames{ind}),'FontSize',7,'FontName','helvetica')
end
xlim([-1 5])
ylim([0 6])
axis off
 saveName = strcat('..\panels\fig_2\fig_S2_ashape_Legend');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
  tprint(saveName,'-SHR -tif',[4 4]);
 %% 
 
load('../data/fig_2/bf_LV.mat')
load('../data/Colours.mat')

colDilution=[1 , 0.8, 0.6 , 0.4];
sz=6;
figure
hold on
nSamp=size(Len_vol_rat,2);

zVec=zeros(nSamp,4);
zVec(:,1)=4*(0:(nSamp-1))+1;
zVec(:,2)=4*(0:(nSamp-1))+2;
zVec(:,3)=4*(0:(nSamp-1))+3;
zVec(:,4)=4*(0:(nSamp-1))+4;


    for bf_ind=1:4
        scatter3(Len_vol_rat(bf_ind,:,1),Len_vol_rat(bf_ind,:,2),zVec(:,bf_ind),sz,colDilution(bf_ind)*Colours(2,:)/256,'filled');
    end

xlim([0 , 8000]);
ylim([0 , 10^7]);
set(gca,'ActivePositionProperty','position','XTick',0:4000:8000,'YTick',[0,5,10]*10^6,'XTickLabel',0:4:8,'YTickLabel',[0,5,10],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');

saveName = strcat('..\panels\fig_2\fig_S2_bbfs');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
  tprint(saveName,'-SHR -tif',[4 4]);

 ShortNames={'bf=0.2','bf=0.4','bf=0.6','bf=0.8'};
Locs=[(2:-1:1)' , zeros(2,1) ; (2:-1:1)' , 3*ones(2,1)];
figure
hold on
for bf_ind=1:4
    scatter(Locs(bf_ind,2),Locs(bf_ind,1),8,colDilution(bf_ind)*Colours(2,:)/256,'filled');
end

for ind=1:4
    text(Locs(ind,2)+0.35,Locs(ind,1),strcat(ShortNames{ind}),'FontSize',7,'FontName','helvetica')
end
xlim([-1 5])
ylim([0 6])
axis off
 saveName = strcat('..\panels\fig_2\fig_S2_bbf_Legend');   
 set(gcf,'renderer','painter','PaperPositionMode','manual');
 
 tprint(saveName,'-SHR -jpg',[4 4]);
 tprint(saveName,'-SHR -eps',[4 4]);
 tprint(saveName,'-SHR -tif',[4 4]);
 
 %% 
load('../data/fig_2/sep_Corrs.mat')
load('../data/Colours.mat')

sz=6;
figure
hold on

colDilutionVec=1-soma_vec(:,1)/100;
nRep=size(soma_vec,1);
colVec=zeros(nRep,3);
for col_ind=1:nRep
    colVec(col_ind,:)=colDilutionVec(col_ind)*Colours(3,:)/256;
end
scatter(soma_vec(:,1),soma_vec(:,2),sz,colVec,'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
[xtrend,ytrend]=scatter_collate(soma_vec(:,1),soma_vec(:,2));
ytrend=smooth(ytrend);

z = zeros(size(xtrend));
col=zeros(2,length(xtrend),3);
for pos_ind=1:length(xtrend)
    col(1,pos_ind,:)=(1-xtrend(pos_ind)/100)*Colours(3,:)/256;
    col(2,pos_ind,:)=(1-xtrend(pos_ind)/100)*Colours(3,:)/256;
end
surface([xtrend;xtrend],[ytrend';ytrend'],[z;z],col,'facecol','no','edgecol','interp','linew',2);

xlim([0 , 100]);
ylim([0 , 6*10^6]);
set(gca,'ActivePositionProperty','position','XTick',0:50:100,'YTick',[0,3,6]*10^6,'XTickLabel',0:50:100,'YTickLabel',[0,3,6],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
%legend('Dendrite','Axon','boxoff','FontSize',6,'FontName','helvetica','Location','northwest');
saveName = strcat('..\panels\fig_2\fig_S2_csepVols');
set(gcf,'renderer','painter','PaperPositionMode','manual');
tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);

figure
hold on
scatter(soma_vec(:,1),soma_vec(:,3),sz,colVec,'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
[xtrend,ytrend]=scatter_collate(soma_vec(:,1),soma_vec(:,3));
ytrend=smooth(ytrend);

z = zeros(size(xtrend));
col=zeros(2,length(xtrend),3);
for pos_ind=1:length(xtrend)
    col(1,pos_ind,:)=(1-xtrend(pos_ind)/100)*Colours(3,:)/256;
    col(2,pos_ind,:)=(1-xtrend(pos_ind)/100)*Colours(3,:)/256;
end
surface([xtrend;xtrend],[ytrend';ytrend'],[z;z],col,'facecol','no','edgecol','interp','linew',2);

xlim([0 , 100]);
ylim([0 , 6000]);
set(gca,'ActivePositionProperty','position','XTick',0:50:100,'YTick',[0,3000,6000],'XTickLabel',0:50:100,'YTickLabel',[0,3,6],'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
%legend('Dendrite','Axon','boxoff','FontSize',6,'FontName','helvetica','Location','northwest');
saveName = strcat('..\panels\fig_2\fig_S2_csepLens');
set(gcf,'renderer','painter','PaperPositionMode','manual');

tprint(saveName,'-SHR -jpg',[4 4]);
tprint(saveName,'-SHR -eps',[4 4]);
tprint(saveName,'-SHR -tif',[4 4]);