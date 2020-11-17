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

soma_seps=[ 10 , 20 , 30 , 40];
nsep=length(soma_seps);

soma_vec=zeros(nsep,nSamp,2);
for soma_ind=1:nsep
    parfor samp_ind=1:nSamp
        test=0;
        while test==0
            try
                nPden=10+round(140*rand(1));
                nPax=10+round(140*rand(1));
                
                
                [dentree] = cube_tree(0.2,nPden,0);
                [axtree] = cube_tree(0.7,nPax,soma_seps(soma_ind));
                
                axtree=resample_tree(axtree,1);
                dentree=resample_tree(dentree,1);
                
                
                [sharedVol,Len1,Len2] = share_boundary_tree(axtree, dentree,'r');
                csyns = peters_tree (axtree, dentree, 2.5,3,'r');
                inSyn=size(csyns,1);
                test=1;
            catch
            end
        end
        soma_vec(soma_ind,samp_ind,:)=[pi/2*Len1*Len2*2.5/sharedVol inSyn];
       [3 soma_ind samp_ind]
    end
end

save('../data/fig_2/sep_syns.mat','soma_vec')
clear