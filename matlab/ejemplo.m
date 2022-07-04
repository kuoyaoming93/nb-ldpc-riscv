%% Borro todas las variables auxiliares
clc, clear;

%% Generacion de matrices H y G
% Cantidad de iteraciones
it_max=10;
% Genero una matriz nula de 837 x it_max
Hdecoded=zeros(it_max,837);
% LDPC matrix (H)
W_base_2
% LDPC Generator Matrix (G)
[G,swap]=NB_LDPC_generator_matrix(gf(H,5)); 

%% Configuracion de parametros iniciales
% LDPC connections for the forward-backward algorithm
[wires]=switches(); 
for i=1:32
    wires_2(i,:)=wires(2,:,i);
end

EbNodB = 6.8:0.1:17;%3.5:0.1:4.5;     % Vector de relacion senial a ruido
seed = [1977 16];                     %[channel data] semilla inicial
LPerrors = 25;                        % Numero de errores a encontrar
Mchannels = 1e10 ;                    % Numero de canales maximo 
rate = 723/837;                       % Code rate
Nbpb = 837*5;                         % Nro de bits por paquete(no coded) 

%% Inicializacion de mas variables y objetos
SNRdB = 10*log10(2*10.^(EbNodB/10)*rate);
tSNR = 'EBN';
tSNRdB = EbNodB;

%% Strings y path para resultados de simulacion
partPrefix = '2(837,723)';   % string with simulation name QCNB-LDPC(n,k)
resultspath='./Results_QCNB_LDPC_in/';
partFile='partFile';

