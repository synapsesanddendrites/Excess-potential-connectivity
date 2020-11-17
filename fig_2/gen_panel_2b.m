% Generates sample morphologies for schematics

nPden=100;

% Generate tree
bf1_tree= cube_tree(0.2,200,0);
bf1_tree= resample_tree(bf1_tree,1);
bf1_tree = jitter_tree(bf1_tree,0.25,25);
bf1_tree = quaddiameter_tree (bf1_tree,0.1); 
bf1_tree.D = bf1_tree.D+1;
bf1_tree=soma_tree(bf1_tree,20);

bf2_tree= cube_tree(0.8,200,0);
bf2_tree= resample_tree(bf2_tree,1);
bf2_tree = jitter_tree(bf2_tree,0.25,25);
bf2_tree = quaddiameter_tree (bf2_tree,0.1); 
bf2_tree.D = bf2_tree.D+1;
bf2_tree=soma_tree(bf2_tree,20);


save('../data/fig_2/morphs_fig_2b.mat')