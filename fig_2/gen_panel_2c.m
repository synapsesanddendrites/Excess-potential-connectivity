% Generates sample morphologies for schematics

nPax=60;
nPden=100;

% Generate trees
dend_tree= cube_tree(0.2,nPden,0);
axonal_tree= cube_tree(0.7,nPax,75);

% Make realistic
dend_tree= resample_tree(dend_tree,1);
axonal_tree=resample_tree(axonal_tree,1);

dend_tree = jitter_tree(dend_tree,0.25,25);
axonal_tree = jitter_tree(axonal_tree,0.25,25);


dend_tree = quaddiameter_tree (dend_tree,0.1); 
axonal_tree = quaddiameter_tree (axonal_tree,0.05); 

dend_tree=cap_tree(dend_tree);
axonal_tree=cap_tree(axonal_tree);

% dend_tree=soma_tree(dend_tree,20);
% axonal_tree=soma_tree(axonal_tree,15);

save('../data/fig_2/morphs_fig_2c.mat','dend_tree','axonal_tree')
clear