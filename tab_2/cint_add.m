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

Trees=CortCol.Trees;
nGroups=length(Groups);

load('../data/fig_3/confdata.mat','paramVec');
load('../data/fig_4/condata.mat')

Permutation=[16,15,6,7,3,2,5,1,4,8,12,13,9,11,10,14];
scats=[];
for axInd=1:size(predCon,1)
    for denInd=1:size(predCon,2)
        nex=length(predCon{axInd,denInd});
            for iex=1:nex
               scats=[scats ; measCon{axInd,denInd}(iex) , pi/2*predCon{axInd,denInd}(iex) , denInd , axInd];
            end
    end
end
scats(isnan(scats))=0;



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
       % try
            iTree=Trees{ind_Class}{ind_Cell};
            iR=iTree.rnames;
            
            Struct=Tab2{ind_Class}{ind_Cell};
            AffVec=Struct.AffVec;
            
            np=length(AffVec);
            
            cInts=[];
            for i_pair=1:np
                in=AffVec(i_pair);
                I=scats(:,1)==in;
                pos_class=scats(I,3)
                
                nZ=length(pos_class);
                if nZ>0
                    fin_scats=[];
                    for c_ind=1:nZ
                        ic=pos_class(c_ind);
                        if ic==Permutation(ind_Class)
                            vals=scats(find(I,c_ind),1);
                            fin_scats=[fin_scats ; vals(end)];
                        end
                    end
                    fin_scats
                    n_poss=length(fin_scats(1,:));
                    n_ind=randi(n_poss);
                    
                    N=fin_scats(2,n_ind);
                    cInt=conf_int(N,0.95,paramVec)
                    cInts=[cInts ; cInt];
                end
            end
            
            Struct.cInts=cInts;
            Tab2{ind_Class}{ind_Cell}=Struct;
            [ind_Class ind_Cell/nCell];
            % catch
            % end
    end
       
end
% save('../data/tab_2/alldata.mat','GmVec','RaVec','Tab2','measSpat') 
% clear