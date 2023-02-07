clc, clear;

nbits=5;
alpha = gf(2,nbits); 

alpha_tb=gf(zeros(1, 2^nbits-1), nbits); 

for i=1:2^nbits-1       
    alpha_tb(i)=alpha^(i-1); 
end

%% GF mult table
gfmult = alpha_tb'*alpha_tb;

%% GF add table
gfadd = repmat([0 alpha_tb],32,1) + repmat([0 alpha_tb],32,1)';

