%SCRIPT PARA GENERAR GRAFICAS AUTOMATICAMENTE

clear all

%PARAMETROS DE LA MATRIZ
N = 32;
M = 16;
K = 16;
q = 16;
p = log2(q);

EBNO_TMM = [0 0.5 1 1.5 2 2.5 3 3.5 4 4.5];
PER_TMM = [1.132813e-001  7.334358e-002 4.107288e-002 1.754464e-002 6.823285e-003 2.243238e-003 4.783769e-004 6.164817e-005 1.187774e-005 1.480310e-006];

%% ET-MM
jj=0;
iter = 50;
clear EC PT

seed = 1977;
EBNO = 0:0.5:4.5;
Nerrors = 25;
ite_max = 50;

Z = ['[C] T-MM LAYERED (' num2str(iter) ' iter), SF=0.9 '];
Y1 = zeros(1,length(EBNO));
for j=1:length(EBNO)
    for i=1:length(seed)
        for k = 1:length(Nerrors)
            archivo2 = ['../RESULTS_ETMM_q16_dc4/(' num2str(N) ',' num2str(K) ')_L' num2str(Nerrors(k)) '[' num2str(seed(i)) ']_N' num2str(N*p) '_K' num2str(K*p) '_EBN' num2str(EBNO(j),'%2.2f') '.txt'];
            if exist(archivo2,'file')
                fid = fopen(archivo2, 'r');
                MNPE_Hdecoded = zeros(1,ite_max);
                MNBE_Hdecoded = zeros(1,ite_max);
                eC = 0;
                nombre = fscanf(fid,'%s', 1);
                for x=1:ite_max
                    MNPE_Hdecoded(x) = fscanf(fid,'%d',1);
                end
                nombre = fscanf(fid,'%s',1);
                for x=1:ite_max
                    MNBE_Hdecoded(x) = fscanf(fid,'%d',1);
                end
                nombre = fscanf(fid,'%s',1);
                eC = fscanf(fid,'%d',1);
                fclose(fid);
                
                Y1(j) = (((MNBE_Hdecoded(iter)/eC)/N)/p);
                break;
            end
        end
    end
end

%% GRAFICA

figure(1), clf
semilogy(EBNO_TMM,PER_TMM,'k--');
hold on;

semilogy(EBNO_TMM,[595/3840 435/4864 418/10624 414/18432 378/51840 324/99840 295/406144 288/2500864 179/14007936 107/(484291*32*4)],'b-*');
semilogy(EBNO,Y1,'r-*');
legend( 'Golden',...
        '[MATLAB] MM (50 iter), SF=1.2 ',...
        Z);
xlabel('E_b/N_o (dB)'), ylabel('FER');
grid on;



