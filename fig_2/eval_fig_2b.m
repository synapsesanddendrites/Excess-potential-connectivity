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

%------------ Balancing factors ------------------------

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
                
                
                [dentree] = cube_tree(balance_facs(bf_ind),nPden,0);
                [axtree] = cube_tree(0.7,nPax,50*rand(1));
                
                axtree=resample_tree(axtree,1);
                dentree=resample_tree(dentree,1);
                
                [sharedVol,Len1,Len2] = share_boundary_tree(axtree, dentree,'r');
                csyns = peters_tree (axtree, dentree, 2.5,3,'r');
                inSyn=size(csyns,1);
                test=1;
            catch
            end
            
        end
        bf_vec(bf_ind,samp_ind,:)=[pi/2*Len1*Len2*2.5/sharedVol inSyn];
        [2 bf_ind samp_ind]
    end
end


%% SIP part
res=50;
ShollVecs=zeros(50,4);
nSampS=1000;
dd=linspace(0.5,400,res);
balance_facs=[ 0.2 , 0.4 , 0.6 , 0.8];
nbf=length(balance_facs);
for bf_ind=1:nbf
    for samp_ind=1:nSampS
        test=0;
        while test==0
           try
                nPden=100;
 
                           
                [dentree] = cube_tree(balance_facs(bf_ind),nPden,0);
                
                [s] = sholl_tree (dentree, dd)';
                ShollVecs(:,bf_ind)=ShollVecs(:,bf_ind)+s/nSampS;
             
                test=1;
            catch
          end
        end
       [2.5 bf_ind samp_ind]
    end
end
save('../data/fig_2/bf_syns.mat','bf_vec','ShollVecs')
clear