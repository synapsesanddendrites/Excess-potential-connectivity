if exist('paramSet','var')==1
    if paramSet==1
    else
        addpath('../functions/trees')
        addpath('../functions')
        start_trees
    end
else
    addpath('../functions/trees')
    addpath('../functions')
    start_trees
end
load('../data/fig_4/cortical_morphs.mat')

%------------ Cortical lengths and volumes ------------------------
Trees=CortCol.Trees;

nClasses=length(Trees);
DendCell=cell(nClasses,1);
AxCell=cell(nClasses,1);

for ind_Class=1:nClasses
    Itrees=Trees{ind_Class};
    nTree=length(Itrees);
    
    
    for TreeInd=1:nTree
        iTree=Itrees{TreeInd};
        
        iR=iTree.rnames;
        for indR=1:length(iR)
            if contains(iR{indR},'Axon')
                axNum=indR;
            end
        end
        axVals=iTree.R==axNum;
        otherVals=iTree.R~=axNum;
        
        iDend=delete_tree (iTree,find(axVals));
        iAxon=delete_tree(iTree,find(otherVals));
        
        ilen=len_tree(iDend);
        iL=sum(ilen(:));
        
        [bound] = boundary_tree(iDend);
        iV=bound.V;
        
        DendCell{ind_Class}(TreeInd,:)=[iL , iV];
        
         ilen=len_tree(iAxon);
        iL=sum(ilen(:));
        
        [bound] = boundary_tree(iAxon);
        iV=bound.V;
        
        AxCell{ind_Class}(TreeInd,:)=[iL , iV];
    end   
end

%% 
Trees=CortCol.Trees;
thresVec=linspace(0.1,0.9,25);

load('../data/fig_4/compdata.mat','GmVec','RaVec')
AllMeans=zeros(16,25);
AllSEMs=zeros(16,25);

for thres_ind=1:25
    thres=thresVec(thres_ind);
   for type_ind=1:16
       nCell=length(Trees{type_ind});
       theseVals=zeros(nCell,1);
       for cell_ind=1:nCell
           iTree=Trees{type_ind}{cell_ind};
           iR=iTree.rnames;
        
        for indR=1:length(iR)
            if contains(iR{indR},'Axon')
                axNum=indR;
            end
        end
        axVals=iTree.R==axNum;
        otherVals=iTree.R~=axNum;
        
        iDend=delete_tree (iTree,find(axVals));
        iDend.Gm=GmVec(type_ind);
        iDend.Ri=RaVec(type_ind);
        
           [Matten]=M_atten_tree(iDend,thres);
           theseVals(cell_ind)=Matten;
       end
       iM=mean(theseVals(:));
       iSE=sqrt(var(theseVals(:))/nCell);
       AllMeans(type_ind,thres_ind)=iM;
       AllSEMs(type_ind,thres_ind)=iSE;
   end
    thres_ind
end
save('../data/fig_3/supdata.mat','thresVec','AllMeans','AllSEMs','AxCell','DendCell')
clear
