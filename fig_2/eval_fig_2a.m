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

%------------ Shapes ------------------------

shape_vec=zeros(4,nSamp,2);
for shape_ind=1:4
    parfor samp_ind=1:nSamp
        test=0;
        while test==0
            try
                nPden=10+round(140*rand(1));
                nPax=10+round(140*rand(1));
                
                if shape_ind==1
                    [dentree] = cube_tree(0.2,nPden,0);
                    [axtree] = cube_tree(0.7,nPax,50*rand(1));
                elseif shape_ind==2
                    [dentree] = sphere_tree(0.2,nPden,200/((4*pi/3)^1/3),0);
                    [axtree] = sphere_tree(0.7,nPax,200/((4*pi/3)^1/3),50*rand(1));
                elseif shape_ind==3
                    [dentree] = cylinder_tree(0.2,nPden,200/((2*pi)^1/3),400/((2*pi)^1/3),0);
                    [axtree] = cylinder_tree(0.7,nPax,200/((2*pi)^1/3),400/((2*pi)^1/3),50*rand(1));
                elseif shape_ind==4
                    [dentree] = cone_tree(0.2,nPden,200*((2*pi/3)^1/3),400*((2*pi/3)^1/3),0);
                    [axtree] = cone_tree(0.7,nPax,200*((2*pi/3)^1/3),400*((2*pi/3)^1/3),50*rand(1));
                end                
                
                axtree=resample_tree(axtree,1);
                dentree=resample_tree(dentree,1);
                
                [sharedVol,Len1,Len2] = share_boundary_tree(axtree, dentree,'r');
                csyns = peters_tree (axtree, dentree, 2.5,3,'r');
                inSyn=size(csyns,1);
                test=1;
                
            catch
            end
        end
        shape_vec(shape_ind,samp_ind,:)=[pi/2*Len1*Len2*2.5/sharedVol inSyn];
        
        [1 shape_ind samp_ind]
    end
end
save('../data/fig_2/shape_syns.mat','shape_vec')
clear