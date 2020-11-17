% Generates sample morphologies for schematics

nPden=200;

% Generate tree
sphere_tree= sphere_tree(0.2,nPden,200/((4*pi/3)^1/3),0);
sphere_tree= resample_tree(sphere_tree,1);
sphere_tree = jitter_tree(sphere_tree,0.25,25);
sphere_tree = quaddiameter_tree (sphere_tree,0.1); 
sphere_tree.D = sphere_tree.D+1;
sphere_tree=soma_tree(sphere_tree,20);

cyl_tree=cylinder_tree(0.2,nPden,200/((2*pi)^1/3),400/((2*pi)^1/3),0);
cyl_tree= resample_tree(cyl_tree,1);
cyl_tree = jitter_tree(cyl_tree,0.25,25);
cyl_tree = quaddiameter_tree (cyl_tree,0.1);
cyl_tree.D = cyl_tree.D+1;
cyl_tree=soma_tree(cyl_tree,20);

cone_tree=cone_tree(0.2,nPden,200*((2*pi/3)^1/3),400*((2*pi/3)^1/3),0);
cone_tree= resample_tree(cone_tree,1);
cone_tree = jitter_tree(cone_tree,0.25,25);
cone_tree = quaddiameter_tree (cone_tree,0.1); 
cone_tree.D = cone_tree.D+1;
cone_tree=soma_tree(cone_tree,20);

save('../data/fig_2/morphs_fig_2a.mat')