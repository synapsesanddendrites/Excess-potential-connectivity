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
load('../data/fig_4/cortical_Rin.mat')

Trees=CortCol.Trees;
nGroups=length(Groups);
spinedis=3;
nClasses=length(Trees);

GmVec=zeros(nClasses,1); %  Get membrane conductivites ---------------
RaVec=zeros(nClasses,1);
for ind_Class=1:nClasses
   [iGm , iRa]=get_Gm_stack(Trees{ind_Class},Rin(ind_Class,1),Rin(ind_Class,2));
    GmVec(ind_Class)=iGm;
    RaVec(ind_Class)=iRa;
    ind_Class
end
%% 

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


Tab2=cell(nClasses,1);

for ind_Class=1:nClasses
    nCell=length(Trees{ind_Class});
    Tab2{ind_Class}=cell(nCell,1);
    
    for ind_Cell=1:nCell
        iTree=Trees{ind_Class}{ind_Cell};
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
        
        NsynAff=0;      
        NsynEff=0;
        NPair=0;
        AffVec=[];
        EffVec=[];
        MpropVec=[];
        
       anyMulti=0;
       whichSlice=[];
        for indmultiGroup=1:lSets
            indGroup=multiGroup(indmultiGroup);
            
            treeInds=Groups{indGroup};
            nsliceTree=size(treeInds,1);
           
         
            for slice_tree_ind=1:nsliceTree
                if treeInds(slice_tree_ind,1)==ind_Class && treeInds(slice_tree_ind,2)==ind_Cell
                    anyMulti=1;
                    whichSlice=indmultiGroup;
                end                
            end
        end
        
        if anyMulti==1
            
            indGroup=multiGroup(whichSlice);
            treeInds=Groups{indGroup};
            nsliceTree=size(treeInds,1);
            
            NPair=nsliceTree-1;       
                  
            for TreeInd=1:nsliceTree
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
                DenTrees{TreeInd}=iDend;
            end
    
            for axInd=1:nsliceTree
                if axInd~=ind_Cell
                    axtree=AxTrees{axInd};
                    csyns = peters_tree (axtree, iDend, spinedis);
                    try
                       [sharedVol,Len1,Len2] = share_boundary_tree(axtree, iDend);
                    catch
                        lTree=len_tree(iDend);
                        Len2=sum(lTree(:));
                    end
                    nsyn=size(csyns,1);
                    NsynAff=NsynAff+nsyn;
                    AffVec=[AffVec ; nsyn];
                    MpropVec=[MpropVec ; Len2];
                end
            end
           
            for denInd=1:nsliceTree
                if denInd~=ind_Cell
                    dentree=DenTrees{denInd};
                    csyns = peters_tree (iAxon, dentree, spinedis);
                    nsyn=size(csyns,1);
                    NsynEff=NsynEff+nsyn;
                    EffVec=[EffVec ; nsyn];
                end
            end
        end
          
        indBP=B_tree(iDend);
        MBranch=nnz(indBP)-1;
        
        iDend.Gm=GmVec(ind_Class);
        iDend.Ri=RaVec(ind_Class);
        MAtten=M_atten_tree(iDend);
        
        lTree=len_tree(iDend);
        Ltree=sum(lTree(:));
        
        Struct.NPair=NPair;
        Struct.NsynAff=NsynAff;
        Struct.NsynEff=NsynEff;
        Struct.AffVec=AffVec;
        Struct.EffVec=EffVec;
        Struct.MpropVec=MpropVec/Ltree;
        Struct.MBranch=MBranch;
        Struct.MAtten=MAtten;  
        
        Tab2{ind_Class}{ind_Cell}=Struct;
        [ind_Class ind_Cell/nCell]
    end

       
end
save('../data/fig_4/compdata.mat','GmVec','RaVec','Tab2','Colours_Cort')
clear