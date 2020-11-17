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

%------------ Length and volume for balancing factors ------------------------

balance_facs=[ 0.2 , 0.4 , 0.6 , 0.8];
nbf=length(balance_facs);

bf_vec=zeros(nbf,nSamp,2);
for bf_ind=1:nbf
    parfor samp_ind=1:nSamp
        test=0;
        while test==0
            try
                nPden=10+round(140*rand(1));
                nPax=10+round(140*rand(1));
                
                
                [dentree] = alex_cube_tree(balance_facs(bf_ind),nPden,0);
                [axtree] = alex_cube_tree(0.7,nPax,50*rand(1));
                
                axtree=resample_tree(axtree,1);
                dentree=resample_tree(dentree,1);
                
                [sharedVol,Len1,Len2] = share_boundary_tree(axtree, dentree,'r');
                csyns = peters_tree (axtree, dentree, 2.5,3,'r');
                inSyn=size(csyns,1);
                test=1;
            catch
            end
            
        end
           
            [bound] = boundary_tree(dentree);
            ilengths=len_tree(dentree);
            Len=sum(ilengths(:));
            
         Len_vol_rat(bf_ind,samp_ind,:)=[Len  bound.V];
        [6 bf_ind samp_ind]
    end
end

save('../data/fig_2/bf_LV.mat','Len_vol_rat')
clear
