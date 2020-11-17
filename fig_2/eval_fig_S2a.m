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

%------------ Length and volume for shapes ------------------------

Len_vol_rat=zeros(4,nSamp,2);
for shape_ind=1:4
    parfor samp_ind=1:nSamp
        test=0;
        while test==0
            try
                nPden=10+round(140*rand(1));
                nPax=10+round(140*rand(1));
                
                if shape_ind==1
                    [dentree] = alex_cube_tree(0.2,nPden,0);
                elseif shape_ind==2
                    [dentree] = alex_sphere_tree(0.2,nPden,200/((4*pi/3)^1/3),0);
                elseif shape_ind==3
                    [dentree] = alex_cylinder_tree(0.2,nPden,200/((2*pi)^1/3),400/((2*pi)^1/3),0);
                elseif shape_ind==4
                    [dentree] = alex_cone_tree(0.2,nPden,200*((2*pi/3)^1/3),400*((2*pi/3)^1/3),0);
                end
             
            [bound] = boundary_tree(dentree);
            ilengths=len_tree(dentree);
            Len=sum(ilengths(:));
                               
           test=1;
           catch
            end
        end
        Len_vol_rat(shape_ind,samp_ind,:)=[Len  bound.V];
       
        
       [5 shape_ind samp_ind]
    end
end
save('../data/fig_2/Shape_LV.mat','Len_vol_rat')
clear
