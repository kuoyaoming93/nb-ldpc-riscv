it_max=10;
Hdecoded=zeros(it_max,837);
%%LDPC matrix
W_base_2
%%LDPC Generator Matrix
[G,swap]=NB_LDPC_generator_matrix(gf(H,5)); 
%%LDPC connections for the forward-backward algorithm
[wires]=switches(); 
for i=1:32
    wires_2(i,:)=wires(2,:,i);
end
EbNodB = 6.8:0.1:17;%3.5:0.1:4.5;     %Signal to noise relation vector
seed = [1977 16];       %[channel data] initial seed
LPerrors = 25;   %number of errors to find 
Mchannels = 1e10 ;  %max number of channels 
rate = 723/837;       %code rate
Nbpb = 837*5;       %number of bits per packet (no coded) 
partPrefix = '2(837,723)'; %string with simulation name QCNB-LDPC(n,k)

% Initialize variables and objects
SNRdB = 10*log10(2*10.^(EbNodB/10)*rate);
tSNR = 'EBN';
tSNRdB = EbNodB;

resultspath='./Results_QCNB_LDPC_in/';
partFile='partFile';

for eS = 1:length(SNRdB)
    tic;
    partFile = [partPrefix '_L' num2str(LPerrors) '[' num2str(seed(1)) '-' num2str(seed(2)) ']_N' num2str(Nbpb) '_K' num2str(Nbpb*rate) '_' tSNR num2str(tSNRdB(eS),'%2.2f') '.mat'];
	try
		clear NPT;
        load([resultspath partFile '.mat']);
        if (exist('NPT','var'))
            display([partFile ' : NPE/NPT = 25/' num2str(MNPT(1,it_max)) ', NBE/NBT = ' num2str(MNBE(1,it_max)) '/' num2str(MNBT(1,it_max)*Nbpb)]);
            %toc;
            continue;
        end
    catch ME
		MNBE_Hdecoded = zeros(1,it_max);
		MNPE_Hdecoded = zeros(1,it_max);
        MNBE = zeros(1,it_max);
        MNPT = zeros(1,it_max);
        MNFL = zeros(1,it_max);
		eC = 1;
	end
    while eC < Mchannels 
        eC;
% 		% Generate and encode a random binary message
        %msg = int32(randint(5,Nbpb*rate/5,2,eC+seed(2)));
        msg = int32(randi(2,5,Nbpb*rate/5)) -1;
        msg_b = 2.^(4:-1:0) * double(msg);

        c = gf(msg_b,5)*G;
        c(swap)=c;
        c_aux=c.x;
        % To check if G and H is correct: c*gf(H',5), the result have to be
        % a 124 zero vector
        % Modulate BPSK
        c = dec2bin(c.x)*1 - 48; %encode bits    
		m = 1 - 2*c; 
        % Noise parameters
        sigma = sqrt(10^(-SNRdB(eS)/10)); 
        %Transmit signal through AWGN channel
        randn('state',eC+seed(1));
        r_ch = m + randn(size(m))*sigma; 
        %r_ch(r_ch<0)=-1; r_ch(r_ch>0)=1;
        %LLR Min_Max
        %reliability_matrx: distance matrix
        %reliability_matrx2: Ln(a)=ln(Pr(cn=sn|channel)/Pr(cn=a|channel))
        
        [reliability_mtrx,reliability_mtrx2]=reliability_pi_7(r_ch,sigma);%%%%%%%%%%%Acelerar
        
        [reliability_mtrx,reliability_mtrx2]=reliability_pi_7(r_ch,sigma);%%%%%%%%%%%Acelerar
        [a,ex]=sort(reliability_mtrx2);
        ex=ex-2; 
        z_HD=gf(2,5).^ex(1,:);
        z_HD(ex(1,:)==-1)=gf(0,5);
        init=length(nonzeros(z_HD~=c_aux))
        %[Q]=L2Q(reliability_mtrx2,H)
        % H: parity check equations (matrix of the LDPC code)
        %Qm,n(a)=Ln(Hm,n*a)
        %%Q for proof and check, comment to obtain the curve
%         [Q]=L2Q(reliability_mtrx2,H);
        %%Good one Q
       
        [Q_2]=L2Q_ver_2(reliability_mtrx2,H); 
        
 %%%%%%%%%%%%%NOT NECESSARY IF FORWARD BACKWARD IS APPLIED%%%%%%%%%%%%%%%%       
        %%Check node processing
        %%Am,n(a):=an'|a+ sum(an')=0
%         [A]=A_group(H,gf_rx,field);    
 %%%%%%%%%%%%%END OF NOT NECESSARY MATLAB FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%% 
 
  for var=1:it_max      
        Rmn=zeros(124*32,837);
        %%Check node processing 
        
        for i=1:124
            [R_Forward,R_Backward,R_Aux,R,col]=R_proof(i,Q_2,H,wires_2);
            Rmn(((i-1)*32+1):(((i-1)*32)+32),col)=R_Aux;
        end
         
        %%Variable node processing
        %%Change Rmn(a) into Rmn(h*a)
        
        Rmn_prime=Rmn_prime_proof(Rmn,H);  
        
        %%Compute: Ln(a)+sum(Rm'n)-min(Qm'n(a))
        
        [Q_prime]=Q_prime_proof(H,Rmn_prime,reliability_mtrx2);
        
        %%Change Q'mn(a) into Q'mn(ha)
        %%BEWARE:ALTHOUGH THE FUNCTION Rmn_prime_proof IS USED WE CHANGE
        %%Q'mn(a) into Q'mn(ha) NOT Rmn
        clear Q_2;
        
        Q_2=Qmn_prime_proof(Q_prime,H); 
        
        %%Tentatively decoding
        %%Compute: Qn(a)=Ln(a)+sum(Rmn)
        
        [Q_n,mag,index]=Qn(H,Rmn_prime,reliability_mtrx2);
        index=index-2;   
        decoded=gf(zeros(1,837),5);
        decoded=gf(2,5).^index;
        decoded(index==-1)=0; 
        
        aux_decoded=decoded.x;
        Hdecoded(var,:) = length(nonzeros(dec2bin(bitxor(uint8(aux_decoded),uint8(c_aux)))*1 - 48));       
        length(nonzeros(aux_decoded~=c_aux))
  end
   

        
%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%         
        %resultado=Hdecoded(1:239)-c2(1:239);
        %The resultado computation is given cause a ZERO CODE is used
%         resultado=index(1:723)-1-c(1:723);
%         if (sum(resultado)~=0) 
%             dknow=1; 
%             error_sym=length(nonzeros((index(1:723)-1)~=c(1:723))); %%%CAUSE ZERO CODE
%         else
%         end
%         %Hdecoded = dec2bin(Hdecoded.x,8)'*1 - 48;    
%         Hdecoded = dec2bin(index(1:723)-1,5)'*1 - 48;   %%%CAUSE ZERO CODE 
        % Calculate BER and PER (partially) 
        %MERR_Hdecoded = sum(abs(double(msg(:))-Hdecoded(:)));
        MERR_Hdecoded = sum(zeros(it_max,837)'+Hdecoded'); %%%CAUSE ZERO CODE 
        MERR_Hdecoded = MERR_Hdecoded/837;   
        MNBE_Hdecoded = MNBE_Hdecoded + MERR_Hdecoded; 
        MNPE_Hdecoded = MNPE_Hdecoded + double(MERR_Hdecoded>0);
        MNBE(MNFL==0) = MNBE_Hdecoded(MNFL==0);
        MNPT(MNFL==0) = eC;
        MNFL = or(MNFL,MNPE_Hdecoded==LPerrors);
        % save partial results
        eC = eC+1; 
        if mod(eC,10)==1
            save([resultspath partFile '.mat'],'MNBE_Hdecoded','MNPE_Hdecoded','MNBE','MNPT','MNFL','eC');
        end
        if (sum(MNFL)==it_max)
            break;
        end
    end
    % save simulation results per SNR
    NPT = eC-1;
    NBT = NPT * Nbpb*rate;
    MNBT = MNPT*Nbpb;
    save([resultspath partFile '.mat'],'MNBE_Hdecoded','MNPE_Hdecoded','MNBE','NBT','MNPT','MNFL','NPT');
    display([partFile ' : NPE/NPT = 25/' num2str(MNPT(1,it_max)) ', NBE/NBT = ' num2str(MNBE(1,it_max)) '/' num2str(MNBT(1,it_max))]);
    toc;
end