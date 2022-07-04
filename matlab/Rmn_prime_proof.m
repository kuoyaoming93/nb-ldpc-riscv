function [Rmn_prime]=Rmn_prime_proof(Rmn,H)
    Rmn_prime=zeros(16*16,32);
%     Rmn_prime_aux=zeros(16*16,32);
%     Raux=zeros(15,1);
%     Raux2=zeros(15,1);
    field=[gf(2,4).^(0:14)];
    field=1./field;
    field=field.x;
    [element,exp]=sort(field);
    exp=double(exp);
    for i=1:16  
        [,col,]=find(H(i,:));
        shift=exp(nonzeros(H(i,:))');
        Baux_2=(0-(15-mod(shift,15)+1));
        Rmn_prime((i-1)*16+1,col)=Rmn((i-1)*16+1,col);
%         Rmn_prime_aux((i-1)*16+1,col)=Rmn((i-1)*16+1,col);
        for b=1: length(col)
            if (Baux_2(b)==1)
                Rmn_prime((i-1)*16+2:(i-1)*16+16,col(b))=[Rmn((i-1)*16+16,col(b)); Rmn((i-1)*16+3-(Baux_2(b)):(i-1)*16+15,col(b))];
            else
                Rmn_prime((i-1)*16+2:(i-1)*16+16,col(b))=[Rmn((i-1)*16+2-(Baux_2(b)):(i-1)*16+16,col(b)); Rmn((i-1)*16+2:(i-1)*16+2-Baux_2(b)-1,col(b)) ];
            end
%             Rmn_prime_aux((i-1)*16+2:(i-1)*16+16,col(b))=circshift(Rmn((i-1)*16+2:(i-1)*16+16,col(b)),[Baux_2(b), 0]);
%             if(nonzeros(Rmn_prime_aux~=Rmn_prime))
%                 stop=1;
%             end
        end
        clear col
    end
end