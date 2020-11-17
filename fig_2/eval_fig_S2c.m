if exist('paramSet','var')==1
    if paramSet==1
    else
        addpath('../functions/trees')
        addpath('../functions')
        start_trees
        nSamp=10000;
    end
else
    addpath('../functions/trees')
    addpath('../functions')
    start_trees
    nSamp=10000;
end


%------------ Intersoma distance ------------------------

soma_vec=zeros(nSamp,3);

parfor samp_ind=1:nSamp
    test=0;
    while test==0
        try
            nPden=10+round(140*rand(1));
            nPax=10+round(140*rand(1));
            
            dist=100*rand(1);
            [dentree] = alex_cube_tree(0.2,nPden,0);
            [axtree] = alex_cube_tree(0.7,nPax,dist);
                        
            [sharedVol,Len1,Len2] = share_boundary_tree(axtree, dentree);
 
            inSyn=size(csyns,1);
            test=1;
        catch
        end
    end
    soma_vec(samp_ind,:)=[dist sharedVol Len1];
    [7 samp_ind]
end


save('../data/fig_2/sep_Corrs.mat','soma_vec')
clear