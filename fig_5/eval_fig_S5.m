addpath('../functions/trees')
addpath('../functions')
start_trees
load('../data/fig_5/raw_morphs.mat')

axNum=2;
nMouse=length(oMouse);
MouseAx=zeros(nMouse,2);
MouseDen=zeros(nMouse,2);
MouseSoma=zeros(nMouse,1);
MouseDens=cell(nMouse,2);

for TreeInd=1:nMouse
    iTree=oMouse{TreeInd};
    
    MouseSoma(TreeInd)=iTree.Y(1);
    
    axVals=iTree.R==axNum;
    otherVals=iTree.R~=axNum;
    
    iDend=delete_tree (iTree,find(axVals));
    iAxon=delete_tree(iTree,find(otherVals));
    
    ilen=len_tree(iDend);
    iL=sum(ilen(:));
    
    [bound] = boundary_tree(iDend);
    iV=bound.V;
    
    MouseDen(TreeInd,:)=[iL , iV];
    
    ilen=len_tree(iAxon);
    iL=sum(ilen(:));
    
    [bound] = boundary_tree(iAxon);
    iV=bound.V;
    
   MouseAx(TreeInd,:)=[iL , iV];
   
   iDendR=resample_tree(iDend,1);
   MouseDens{TreeInd,1}=iDendR.Y;
   iAxR=resample_tree(iAxon,1);
   MouseDens{TreeInd,2}=iAxR.Y;
   
end

nHuman=length(oHuman);
HumanAx=zeros(nHuman,2);
HumanDen=zeros(nHuman,2);
HumanSoma=zeros(nHuman,1);
HumanDens=cell(nHuman,2);
for TreeInd=1:nHuman
    iTree=oHuman{TreeInd};
    
    HumanSoma(TreeInd)=iTree.Y(1);
  
    axVals=iTree.R==axNum;
    otherVals=iTree.R~=axNum;
    
    iDend=delete_tree (iTree,find(axVals));
    iAxon=delete_tree(iTree,find(otherVals));
    
    ilen=len_tree(iDend);
    iL=sum(ilen(:));
    
    [bound] = boundary_tree(iDend);
    iV=bound.V;
    
    HumanDen(TreeInd,:)=[iL , iV];
    
    ilen=len_tree(iAxon);
    iL=sum(ilen(:));
    
    [bound] = boundary_tree(iAxon);
    iV=bound.V;
    
   HumanAx(TreeInd,:)=[iL , iV];
   
    iDendR=resample_tree(iDend,1);
   HumanDens{TreeInd,1}=iDendR.Y;
   iAxR=resample_tree(iAxon,1);
   HumanDens{TreeInd,2}=iAxR.Y;
end

save('../data/fig_5/supdata.mat','MouseAx','MouseDen','MouseDens','MouseSoma','HumanAx','HumanDen','HumanDens','HumanSoma') 
clear
