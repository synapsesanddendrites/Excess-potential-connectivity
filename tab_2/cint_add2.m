addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/fig_4/cortical_morphs.mat')
load('../data/fig_4/cortical_Rin.mat')

Trees=CortCol.Trees;
nGroups=length(Groups);
spinedis=3;
nClasses=length(Trees);

load('../data/tab_2/alldata.mat','GmVec','RaVec','Tab2','measSpat') 
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


for ind_Class=1:nClasses
    nCell=length(Trees{ind_Class});
    
    for ind_Cell=1:nCell
        Struct=Tab2{ind_Class}{ind_Cell}
        AffVec=Struct.AffVec;
        
        np=length(AffVec);
        cInts=[];
        
        if np>0
            
            for i_pair=1:np
                in=AffVec(i_pair);
                cInt=conf_int(in,95,paramVec);
                cInts=[cInts ; cInt];
            end
        end
        
        Struct.cInts=cInts;
        Tab2{ind_Class}{ind_Cell}=Struct;
        
        
        [ind_Class ind_Cell/nCell]
    end
    
    
end