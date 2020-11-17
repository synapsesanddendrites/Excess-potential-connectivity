%% Trend in terms of volume
 clear

nstep=10;
load('../data/fig_1/dendrite_panel.mat')
load('../data/fig_1/axon_panel.mat')
load('../data/fig_1/volume_panel.mat')
load('../data/Colours.mat')

lD=2.4e+03;
lA=3.0e+03;
vS=2.4e+06;
spineVec=[1.5,2.5,3.5,4.5];

names={'denVals','axVals','shVals'};
alltitles={'Dendrite lengths (variable)','Axon lengths (fixed)','Volumes (fixed)';'Dendrite lengths (fixed)','Axon lengths (variable)','Volumes (fixed)';'Dendrite lengths (fixed)','Axon lengths (fixed)','Volumes (variable)';'Dendrite lengths (fixed)','Axon lengths (fixed)','Volumes (fixed)'};
reorder=[2,3,1];

dColours=[Colours(1,:)/256 ; Colours(2,:)/256 ; 0.3 , 0.3 , 0.3];
%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%
figure

for factor=1:3
    iname=names{factor};
    iDists=eval(iname);
    for type=1:3
        theseVals=iDists(:,reorder(type));
        
        minVal=min(theseVals(:));
        maxVal=max(theseVals(:));
        
        edges=linspace(minVal,maxVal,nstep+1);
        theseDens=histcounts(theseVals,edges)/length(theseVals);
        
        delVal=(maxVal-minVal)/(nstep+1);
        mids=linspace(minVal+delVal/2,maxVal-delVal/2,nstep);
        
        subplot(3,3,(factor-1)*3+type)
        hold on
        bar(mids,theseDens,'FaceColor',dColours(type,:),'EdgeColor',dColours(type,:))
        maxDen=max(theseDens(:));
        
        ylim([0 0.4])
        if type==1
            scatter(lD,1.1*maxDen,4,[0,0,0],'filled')
            if factor==1              
                xlim([0 6000])               
                set(gca,'ActivePositionProperty','position','XTick',0:2000:6000,'XTickLabel',0:2:6,'YTick',0:0.2:0.4,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
            else
                 xlim([2200 2600])               
                set(gca,'ActivePositionProperty','position','XTick',2200:200:2600,'XTickLabel',2.2:0.2:2.6,'YTick',0:0.2:0.4,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
            end
            elseif type==2
                scatter(lA,1.1*maxDen,6,[0,0,0],'filled')
                if factor==2
                    xlim([0 6000])
                    set(gca,'ActivePositionProperty','position','XTick',0:2000:6000,'XTickLabel',0:2:6,'YTick',0:0.2:0.4,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
                else
                    xlim([2800 3200])
                    set(gca,'ActivePositionProperty','position','XTick',2800:200:3200,'XTickLabel',2.8:0.2:3.2,'YTick',0:0.2:0.4,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
                end           
        elseif type==3
            scatter(vS,1.1*maxDen,6,[0,0,0],'filled')
             if factor==3
                    xlim([1500000 3500000])
                    set(gca,'ActivePositionProperty','position','XTick',1500000:500000:3000000,'XTickLabel',0:1:3,'YTick',0:0.2:0.4,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
                else
                    xlim([2200000 2600000])
                    set(gca,'ActivePositionProperty','position','XTick',2200000:100000:2600000,'XTickLabel',2.2:0.1:2.6,'YTick',0:0.2:0.4,'ticklength',[0.04 0.08],'XMinorTick','on','YMinorTick','on','fontsize',8,'fontname','helvetica','tickdir','out');
                end           
            
        end
       % title(alltitles{factor,type});
    end
    
end
saveName='../panels/fig_1/fig_S1';
tprint(saveName,'-SHR -jpg',[12 12]);
set(gcf,'renderer','painter');
tprint(saveName,'-SHR -eps',[12 12]);
tprint(saveName,'-SHR -tif',[12 12]);