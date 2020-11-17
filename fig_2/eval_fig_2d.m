if exist('paramSet','var')==1
    if paramSet==1
    else
        addpath('../functions/trees')
        addpath('../functions')
        start_trees
        nSamp=50000;
    end
else
    addpath('../functions/trees')
    addpath('../functions')
    start_trees
    nSamp=50000;
end

%------------ DSCAM null ------------------------


iterations 	= [0 100 400 700];
nit=length(iterations);
DSCAM_vec=zeros(nSamp,nit,4);
Len_vol_rat=zeros(nSamp,nit,2);
parfor samp_ind=1:nSamp
    test=0;
    while test==0
        
        nPden=10+round(140*rand(1));
        nPax=10+round(140*rand(1));
        
        [outtree] = cube_tree(0.2,nPden,0);
        [axtree] = cube_tree(0.7,nPax,100*rand(1));
        
        try
            [sharedVol,Len1,Len2] = share_boundary_tree(axtree, outtree);
            csyns = peters_tree (axtree, outtree, 2.5);
            nSyn=size(csyns,1);
            iDS=zeros(nit,4);
            iLV=zeros(nit,2);
            
             
             [bound] = boundary_tree(outtree);
             ilengths=len_tree(outtree);
             Len=sum(ilengths(:));
             iLV(1,:)=[Len bound.V];
            
              iDS(1,:)=[Len bound.V pi/2*Len1*Len2*2.5/sharedVol nSyn];

            for it_ind=2:length(iterations)
                
                outtree = dscam_tree(outtree, iterations(it_ind));
                
                csyns = peters_tree (axtree, outtree, 2.5);
                nSyn=size(csyns,1);
            
                
                [bound] = boundary_tree(outtree);
             ilengths=len_tree(outtree);
             Len=sum(ilengths(:));
                iLV(it_ind,:)=[Len bound.V];

            iDS(it_ind,:)=[Len bound.V pi/2*Len1*Len2*2.5/sharedVol nSyn];

            end
            DSCAM_vec(samp_ind,:,:)=iDS;
            Len_vol_rat(samp_ind,:,:)=iLV;
            [4 samp_ind]
            test=1;
       catch
        end
    end
end
save('../data/fig_2/DSCAM_syns.mat','Len_vol_rat','DSCAM_vec')
clear