if exist('paramSet','var')==1
    if paramSet==1
    else
        addpath('../functions/trees')
        addpath('../functions')
        start_trees
        nSamp=1000;
        
        lD=2.4e+03;
        lA=3.0e+03;
        vS=2.4e+06;
        spineVec=[1.5,2.5,3.5,4.5];
    end
else
    addpath('../functions/trees')
    addpath('../functions')
    start_trees
    nSamp=1000;
    
    lD=2.4e+03;
    lA=3.0e+03;
    vS=2.4e+06;
    spineVec=[1.5,2.5,3.5,4.5];
end


%------------ Axon Unknown ------------------------
axVals=zeros(nSamp,3);
axTrend=zeros(nSamp,2,length(spineVec));
parfor iRep=1:nSamp
    test=0;
    while test==0
        try
            nAll=0;
            [outstruct] = genTree_La(lD,vS);
            axVals(iRep,:)=[outstruct.shVol , outstruct.Ld , outstruct.La];
            
            iTrend=zeros(2,length(spineVec));
            for iSpine=1:length(spineVec)
                csyns = peters_tree (outstruct.dentree, outstruct.axtree, spineVec(iSpine),3,'r');
                nSyn=size(csyns,1);
                iTrend(1,iSpine)=outstruct.La;
                iTrend(2,iSpine)=nSyn;
            end
            axTrend(iRep,:,:)=iTrend;
            if outstruct.error<=0.05
                test=1;
                [1 run_ind iRep]
            end
        catch
        end
    end
end
savename='../data/fig_1/axon_panel.mat';
save(savename,'axTrend','axVals')
clear