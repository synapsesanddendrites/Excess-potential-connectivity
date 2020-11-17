addpath('../functions/trees')
addpath('../functions')
start_trees

load('../data/fig_4/cortical_morphs.mat')
load('../data/fig_4/cortical_Rin.mat')

Trees=CortCol.Trees;
nGroups=length(Groups);
spinedis=3;
nClasses=length(Trees);

%  Get membrane conductivites ---------------
RaVec=zeros(nClasses,1);
GmVec=zeros(nClasses,1); 
for ind_Class=1:nClasses
   [iGm , iRa]=get_Gm_stack(Trees{ind_Class},Rin(ind_Class,1),Rin(ind_Class,2));
    GmVec(ind_Class)=iGm;
    RaVec(ind_Class)=iRa;
    ind_Class
end
%% 
load('../data/fig_3/confdata.mat','paramVec');
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
        cInts=[];
        
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
                    try
                        csyns = peters_tree (axtree, iDend, spinedis);
                    catch
                        csyns=[];
                    end
                    try
                        [sharedVol,Len1,Len2] = share_boundary_tree(axtree, iDend);
                        cInt=conf_int(pi/2*Len1*Len2*3/sharedVol,95,paramVec);
                    catch
                        lTree=len_tree(iDend);
                        Len2=sum(lTree(:));
                        cInt=conf_int(size(csyns,1),0.95,paramVec);
                    end                 
                    nsyn=size(csyns,1);
                    NsynAff=NsynAff+nsyn;
                    AffVec=[AffVec ; nsyn];
                    cInts=[cInts ; cInt];
                    MpropVec=[MpropVec ; Len2];
                end
            end
           
            for denInd=1:nsliceTree
                if denInd~=ind_Cell
                    dentree=DenTrees{denInd};
                    try
                        csyns = peters_tree (iAxon, dentree, spinedis);
                    catch
                        csyns=[];
                    end
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
        Struct.cInts=cInts;
        
        Tab2{ind_Class}{ind_Cell}=Struct;
        [ind_Class ind_Cell/nCell]
    end

       
end
%% 

%%%%%%%%%%%%%%%%%%%%%%%%% Nearest neighbour ratios %%%%%%%%%%%%%%%%%%%%%
spinedis=3;

nClasses=length(Trees);
measSpat=cell(nClasses);

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
    
    %     treeInds=Groups{indGroup}
    %     nTree=size(treeInds,1)
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

            else
                 dentree=DenTrees{denInd};
                 
                 [nnrObj] = synapse_nnr(axtree, dentree, spinedis);
                                
                 
                 theseMeas=measSpat{AxTypes(axInd),DenTypes(denInd)};
                 theseMeas=[theseMeas ; nnrObj.nnr , nnrObj.ext , nnrObj.anr , nnrObj.extA];
                 measSpat{AxTypes(axInd),DenTypes(denInd)}=theseMeas;
                 
            end
            catch
            end
            
        end
    end

       
end
save('../data/tab_2/alldata.mat','GmVec','RaVec','Tab2','measSpat') 
clear