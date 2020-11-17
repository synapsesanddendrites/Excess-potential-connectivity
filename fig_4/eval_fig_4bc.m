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
predCon=cell(nClasses);
measCon=cell(nClasses);
spineDists=[];

multiDays=zeros(nGroups,1); % Days 	multiple morphologies

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
        for denInd=1:nTree
            if denInd==axInd
                % Autapses!
            else
                dentree=DenTrees{denInd};
                csyns = peters_tree (axtree, dentree, spinedis);
                nsyn=size(csyns,1);
                if nsyn>0
                    spineDists=[spineDists ; csyns(:,3)];
                end
                
                theseMeas=measCon{AxTypes(axInd),DenTypes(denInd)};
                theseMeas=[theseMeas ; nsyn];
                measCon{AxTypes(axInd),DenTypes(denInd)}=theseMeas;
                
                try
                    [sharedVol,Len1,Len2] = share_boundary_tree(axtree, dentree);
                    nsynEst=Len1*Len2*spinedis/sharedVol;
                catch
                    nsynEst=nsyn;
                end
                
                thesePreds=predCon{AxTypes(axInd),DenTypes(denInd)};
                 thesePreds=[thesePreds ; nsynEst];
                 predCon{AxTypes(axInd),DenTypes(denInd)}=thesePreds;
                 
            end
            
        end
    end

       
end
save('../data/fig_4/condata.mat','Colours_Cort','predCon','measCon') 
clear