old=cd('./fig_1');
save('Current_directory','old');
eval_fig_1
load('Current_directory','old');
delete('Current_directory');
cd(old)
%%
old=cd('./fig_2');
save('Current_directory','old');
eval_fig_2
load('Current_directory','old');
delete('Current_directory');
cd(old)
%%
old=cd('./fig_3');
save('Current_directory','old');
eval_fig_3
load('Current_directory','old');
delete('Current_directory');
cd(old)
%%
old=cd('./fig_4');
save('Current_directory','old');
eval_fig_4
load('Current_directory','old');
delete('Current_directory');
cd(old)
%%
old=cd('./fig_5');
save('Current_directory','old');
eval_fig_5
load('Current_directory','old');
delete('Current_directory');
cd(old)
%%
old=cd('./tab_2');
save('Current_directory','old');
eval_tab_2
load('Current_directory','old');
delete('Current_directory');
cd(old)
clear