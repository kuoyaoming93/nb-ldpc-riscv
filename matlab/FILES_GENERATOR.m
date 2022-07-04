%% Matrix Shu-Lin Optical
ShuLin=0;

%% Code parameters
dc=4; 
dv=2; 
%q_QC=16;
%qQC_REAL = 21;

%m_QC=log2(q_QC);
%if ShuLin==0
%    N_code=(q_QC)*(dc); 
%    M_code=(q_QC)*dv; 
%else
%    N_code=(q_QC-1)*(dc+1); 
%    M_code=(q_QC-1)*dv; 
%end

N_code = 32;
M_code = 16;

%% Field parameters
q_field=16; 
m_field=log2(q_field);


%% GF TABLES 

%GF add table
gfadd=repmat([0 gf(2,m_field).^(0:(q_field-2))],q_field,1)+repmat([0 gf(2,m_field).^(0:(q_field-2))]',1,q_field);
gfadd=gfadd.x; 
gfadd=uint32(gfadd);
[x, a] = sort(gfadd(1,:));
[x, b] = sort(gfadd(:,1));

gfadd_temp = zeros(q_field);
for i=1:q_field
    gfadd_temp(i,:) = gfadd(i,a);
end
for i=1:q_field
    gfadd_temp(:,i) = gfadd_temp(b,i);
end


%Elements table sort by integer equivalent
fld=gf(2,m_field).^(0:q_field-2); 
fld=fld.x;

%Elements table sort by power equivalent
[ a , pow]=sort(fld);
pow=pow-1;
pow=[0 pow];

%Converts from MATLAB position to exp on GF(q)
pos_exp_table=gfadd(:,1); 

%Convierte de exp+1 en GF(q) a posicion en MATLAB
exp_pos_table=zeros(q_field,1);
exp_pos_table(pos_exp_table+1)=1:q_field; 
exp_pos_table = uint32(exp_pos_table);


%GF table for BPSK modulation
fld_bpsk=gf(zeros(1,(q_field-1)));
fld_bpsk(2:q_field)=gf(2,m_field).^(0:(q_field-2));
fld_bpsk=dec2bin(fld_bpsk.x)*1-48;

fld_bpsk2 = fld_bpsk(exp_pos_table,:);


a = gf( 1:q_field-1, m_field);      
table = a' * a;        
gfmult = zeros(q_field);
gfmult(2:end,2:end) = table.x;

%Creacion del archivo gf_tables.h
dir = '../';
nombre = ['gf_tables' '_GF' num2str(q_field)];
f = sprintf([dir nombre '.h']);
pack = fopen(f,'w');

    fprintf(pack,['/* GENERATED WITH MATLAB...*/' '\n']);
    fprintf(pack,'\n');
    fprintf(pack,['/* Tables for Finite Field GF(' num2str(q_field) ')*/\n']);
    fprintf(pack,'\n');
    
    fprintf(pack,['const int gfadd[' num2str(q_field) '][' num2str(q_field) '] = { {']);
    for i=1:q_field
        for j=1:q_field
            fprintf(pack,num2str(gfadd_temp(i,j)));
            if (j~=q_field)
                fprintf(pack,',');
            else
                fprintf(pack,'}');
                if (i~=q_field)
                    fprintf(pack,', \n                          {');
                else
                    fprintf(pack,'}; \n \n');
                end
            end
        end
    end
    
    fprintf(pack,['const int gfmult[' num2str(q_field) '][' num2str(q_field) '] = { {']);
    for i=1:q_field
        for j=1:q_field
            fprintf(pack,num2str(gfmult(i,j)));
            if (j~=q_field)
                fprintf(pack,',');
            else
                fprintf(pack,'}');
                if (i~=q_field)
                    fprintf(pack,', \n                          {');
                else
                    fprintf(pack,'}; \n \n');
                end
            end
        end
    end
    
    fprintf(pack,['const int fld_bpsk[' num2str(q_field) '][' num2str(m_field) '] = { {']);
    for i=1:q_field
        for j=1:m_field
            fprintf(pack,num2str(fld_bpsk2(i,j)));
            if (j~=m_field)
                fprintf(pack,',');
            else
                fprintf(pack,'}');
                if (i~=q_field)
                    fprintf(pack,', \n                             {');
                else
                    fprintf(pack,'}; \n \n');
                end
            end
        end
    end
    
    fprintf(pack,['const int fld_bpsk_bip[' num2str(q_field) '][' num2str(m_field) '] = { {']);
    for i=1:q_field
        for j=1:m_field
            fprintf(pack,num2str(1-2*fld_bpsk2(i,j)));
            if (j~=m_field)
                fprintf(pack,',');
            else
                fprintf(pack,'}');
                if (i~=q_field)
                    fprintf(pack,', \n                             {');
                else
                    fprintf(pack,'}; \n \n');
                end
            end
        end
    end
    
    fprintf(pack,['const int fld_bpsk_order[' num2str(q_field) '][' num2str(m_field) '] = { {']);
    for i=1:q_field
        for j=1:m_field
            fprintf(pack,num2str(fld_bpsk(i,j)));
            if (j~=m_field)
                fprintf(pack,',');
            else
                fprintf(pack,'}');
                if (i~=q_field)
                    fprintf(pack,', \n                             {');
                else
                    fprintf(pack,'}; \n \n');
                end
            end
        end
    end
    
fclose(pack);


%% CHECK NODE TABLE SPECIFICS FOR EACH ALGORITHM 

    %Converts from MATLAB position to exp on GF(q)
    pos_exp_table=gfadd(:,1); 
    
    %Convierte de exp+1 en GF(q) a posicion en MATLAB
    exp_pos_table=zeros(q_field,1);
    exp_pos_table(pos_exp_table+1)=1:q_field; 
    exp_pos_table = uint32(exp_pos_table);

    %Tabla para las configuraciones con dos desviaciones
    conf_tb=zeros(q_field/2-1,2,q_field-1);
    for i=1:q_field-1
        [row col]=find(gfadd==pos_exp_table(i+1));
        matriz = [row';col']';
        matriz=sort(matriz,2);
        matriz=sortrows(matriz,1);
        matriz=matriz(3:2:end,:);
        conf_tb(:,:,i)=matriz;
    end
    conf_tb=uint32(conf_tb);

    %Creacion del archivo CNU_tables.h
    dir = '../';
    nombre = ['CNU_tables' '_GF' num2str(q_field)];
    f = sprintf([dir nombre '.h']);
    pack = fopen(f,'w');

        fprintf(pack,['/* GENERATED WITH MATLAB...*/' '\n']);
        fprintf(pack,'\n');
        fprintf(pack,['/* CNU tables for Finite Field GF(' num2str(q_field) ')*/\n']);
        fprintf(pack,'\n');

        fprintf(pack,['const int exp_pos_table[' num2str(q_field) '] = {']);
        for j=1:q_field
            fprintf(pack,num2str(exp_pos_table(j)));
            if (j~=q_field)
                fprintf(pack,',');
            else
                fprintf(pack,'}; \n \n');
            end
        end
        
        fprintf(pack,['const int pos_exp_table[' num2str(q_field) '] = {']);
        for j=1:q_field
            fprintf(pack,num2str(pos_exp_table(j)));
            if (j~=q_field)
                fprintf(pack,',');
            else
                fprintf(pack,'}; \n \n');
            end
        end
        
        fprintf(pack,['const int conf_tb[' num2str(q_field/2-1) '][2][' num2str(q_field-1) '] = { { {']);
        for i=1:q_field/2-1
            for j=1:2
                for k=1:q_field-1
                    fprintf(pack,num2str(conf_tb(i,j,k)));
                    if (k~=q_field-1)
                        fprintf(pack,',');
                    else
                        fprintf(pack,'}');
                        if (j~=2)
                            fprintf(pack,', \n                                 {');
                        else
                            fprintf(pack,'}');
                            if (i~=q_field/2-1)
                                fprintf(pack,', \n\n                                 {{');
                            else
                                fprintf(pack,'}; \n \n');
                            end                        
                        end
                    end
                end
            end
        end

    fclose(pack);    
    
    
%% PARITY CHECK MATRIX

    %LDPC parity check matrix
%     W_base_optical_RS_concat2;
    H_matrix;
    %load H_150_600_GF256.mat;
    
%     HH = zeros(dv*qQC_REAL,dc*qQC_REAL);
%     for jj = 1:dv
%       for ii = 1:dc
%          HH(qQC_REAL*(jj-1)+1:qQC_REAL*(jj),qQC_REAL*(ii-1)+1:qQC_REAL*(ii)) = ...
%             H(q_QC*(jj-1)+1:q_QC*(jj-1)+qQC_REAL,q_QC*(ii-1)+1:q_QC*(ii-1)+qQC_REAL);
%       end
%     end
%     H=HH;
%     
%     for ii = 1:dv*qQC_REAL
%        dc_list(ii) = length(find(H(ii,:)));
%     end
%     dc = max(dc_list);
% 
%    for ii = 1:dc*qQC_REAL
%        dv_list(ii) = length(find(H(:,ii)));
%     end
%     dv = max(dv_list);
%     M_code = dv*qQC_REAL;
%     N_code = (dc+1)*qQC_REAL;
    
    %Table generation for coef(H) and non zero position of H matrix
    col = zeros(M_code,dc); % Non-zero element positions for each row of H
    coefH = zeros(size(col)); % Non-zero element power values for each row of H
    for i=1:M_code
%    for i=1:dv*qQC_REAL
        temp = find(H(i,:));
        col(i,1:length(temp))= temp;
        temp = nonzeros(H(i,:));
        coefH(i,1:length(temp))= temp;
    end

    pow_coefH = uint32(pow(coefH+1));

    %LDPC Generator Matrix
    [G,swap]=NB_LDPC_generator_matrix(gf(H,m_field)); 
    GG = G.x;
    K_code_aux=size(G);
    K_code=K_code_aux(1);
    rate = K_code/N_code;       	% code rate
    Nbpb = N_code*m_field;       	% number of bits per packet (no coded) 
    
    %msg = randi([0 1],m_field,Nbpb*rate/m_field,'int32');
    msg = randi([0 1],m_field,K_code,'int32');
    msg_b = 2.^((m_field-1):-1:0) * double(msg);
    c = gf(msg_b,m_field)*G;
    c(swap)=c;
    c_aux=c.x;

    % Modulate BPSK
    c = dec2bin(c.x)*1 - 48; %encode bits    
    m = 1 - 2*c;

    %Creacion del archivo que almacena la informacion de la matriz
    dir = '';
    nombre = ['Hmat_N' num2str(N_code) '_M' num2str(M_code) '_GF' num2str(q_field)];
    f = sprintf([dir nombre '.h']);
    pack = fopen(f,'w');

        fprintf(pack,['/* GENERATED WITH MATLAB...*/' '\n']);
        fprintf(pack,'\n');
        fprintf(pack,['/* Matrix parameters for the code (' num2str(N_code) ',' num2str(M_code) ') over GF(' num2str(q_field) ')*/\n']);
        fprintf(pack,'\n');
        
        fprintf(pack,'/* Code Parameters */ \n');
        fprintf(pack,['#define q_field ' num2str(q_field) '\n']);
        fprintf(pack,['#define m_field ' num2str(m_field) '\n']);
        fprintf(pack,['#define N_code ' num2str(N_code) '\n']);
        fprintf(pack,['#define M_code ' num2str(M_code) '\n']);
        fprintf(pack,['#define dc ' num2str(dc) '\n']);
        fprintf(pack,['#define dv ' num2str(dv) '\n']);
%        fprintf(pack,['#define q_QC ' num2str(q_QC) '\n \n']);
        
        fprintf(pack,['#define K_code ' num2str(K_code) '\n']);
        fprintf(pack,['#define rate ' num2str(rate) '\n']);
        fprintf(pack,['#define Nbpb ' num2str(Nbpb) '\n \n']);
        
        fprintf(pack,'/* Non-zero element positions for each row of H */ \n');
        fprintf(pack,['const int col[' num2str(M_code) '][' num2str(dc) '] = { {']);
        for i=1:M_code
            for j=1:dc
                fprintf(pack,num2str(col(i,j)-1));
                if (j~=dc)
                    fprintf(pack,',');
                else
                    fprintf(pack,'}');
                    if (i~=M_code)
                        fprintf(pack,', \n                          {');
                    else
                        fprintf(pack,'}; \n \n');
                    end
                end
            end
        end

        fprintf(pack,'/* Non-zero element power values for each row of H */ \n');
        fprintf(pack,['const int pow_coefH[' num2str(M_code) '][' num2str(dc) '] = { {']);
        for i=1:M_code
            for j=1:dc
                fprintf(pack,num2str(pow_coefH(i,j)));
                if (j~=dc)
                    fprintf(pack,',');
                else
                    fprintf(pack,'}');
                    if (i~=M_code)
                        fprintf(pack,', \n                          {');
                    else
                        fprintf(pack,'}; \n \n');
                    end
                end
            end
        end
        
        fprintf(pack,['const int m[' num2str(N_code) '][' num2str(m_field) '] = { {']);
        for i=1:N_code
        %for i=1:dc*qQC_REAL
            for j=1:m_field
                fprintf(pack,num2str(m(i,j)));
                if (j~=m_field)
                    fprintf(pack,',');
                else
                    fprintf(pack,'}');
                    if (i~=N_code)
                        fprintf(pack,', \n                             {');
                    else
                        fprintf(pack,'}; \n \n');
                    end
                end
            end
        end

        

    fclose(pack);        
    
    %save H_N504_K252_GF16 col coefH G swap pow_coefH gfadd exp_pos_table pos_exp_table conf_tb
    save(['Hmat_N' num2str(N_code) '_M' num2str(M_code) '_GF' num2str(q_field) '.mat'], 'col', 'coefH', 'G', 'swap', 'pow_coefH', 'gfadd', 'exp_pos_table', 'pos_exp_table', 'conf_tb');
   
    
