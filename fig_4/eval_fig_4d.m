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

Trees=CortCol.Trees;
nGroups=length(Groups);
spinedis=3;

nClasses=length(Trees);
measSpat=cell(nClasses);
vol_rats=cell(nClasses);

multiDays=zeros(nGroups,1); % Days with multiple morphologies

for indGroup=1:nGroups
     treeInds=Groups{indGroup};
    
    nTree=size(treeInds,1);
    if nTree>1
        multiDays(indGroup)=1;
    end
    
end

lSets=nnz(multiDays);
multiGroup=find(multiDays);

for indmultiGroup=1:lSets
 
    indGroup=multiGroup(indmultiGroup);

    treeInds=Groups{indGroup};
    nTree=size(treeInds,1);
    
    
    AxTrees=cell(nTree,1);
    AxTypes=zeros(nTree,1);
    
    DenTrees=cell(nTree,1);
    DenTypes=zeros(nTree,1);
    
    for TreeInd=1:nTree
        iTree=Trees{treeInds(TreeInd,1)}{treeInds(TreeInd,2)};
        
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
        
        AxTrees{TreeInd}=iAxon;
        AxTypes(TreeInd)=treeInds(TreeInd,1);
        
        DenTrees{TreeInd}=iDend;
        DenTypes(TreeInd)=treeInds(TreeInd,1);
    end
    
    for axInd=1:nTree
        axtree=AxTrees{axInd};
        denInd=1;
        for denInd=1:nTree
            try
            if denInd==axInd
                % Autapses!
            else
                 dentree=DenTrees{denInd};
                 
                 [nnrObj] = synapse_nnr(axtree, dentree, spinedis);
                                
                 
                 theseMeas=measSpat{AxTypes(axInd),DenTypes(denInd)};
                 theseMeas=[theseMeas ; nnrObj.nnr , nnrObj.ext , nnrObj.anr , nnrObj.extA];
                 measSpat{AxTypes(axInd),DenTypes(denInd)}=theseMeas;
                 
                 
                 [volObj] = synapse_volr(axtree, dentree, spinedis);                                              
                 theseMeas=vol_rats{AxTypes(axInd),DenTypes(denInd)};
                 theseMeas=[theseMeas ; volObj.volr , volObj.ext];
                 vol_rats{AxTypes(axInd),DenTypes(denInd)}=theseMeas;
                 
            end
            catch
            end
            
        end
    end

       
end
save('../data/fig_4/spatdata.mat','measSpat','vol_rats','Colours_Cort')
clear