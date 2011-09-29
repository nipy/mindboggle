% test find root node
clear;clc;close all;
path(path,'toolbox') ;
% sk_filename='../result/cylinder1_contract_t(2)_nn(14)_WL(10.633697)_WH(1.000000)_sl(3.000000)_skeleton.mat';
% sk_filename='../result/simplejoint_v4770_contract_t(3)_nn(30)_WL(15.378798)_WH(1.000000)_sl(3.000000)_skeleton.mat';
sk_filename='../result/horse_v1987_contract_t(3)_nn(24)_WL(7.786614)_WH(1.000000)_sl(3.000000)_skeleton.mat';

load(sk_filename,'M');

%%
[joints, segments] = find_joints(M, false);
[root_id, global_dist, skel_size] = find_root_node(M, joints, true);