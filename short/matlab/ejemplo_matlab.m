C = 2;          
K = 2560;       
F = 36;    
txcbs = ones(K-F,C);
fillers = -1*ones(F,C);
txcbs = [txcbs;fillers];  

bgn = 2; 
txcodedcbs = nrLDPCEncode(txcbs,bgn);