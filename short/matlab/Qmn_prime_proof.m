function [Qmn_prime]=Qmn_prime_proof(Qmn,H)
    Qmn_prime=zeros(16*16,32);
    field=[gf(2,4).^(0:14)];
    field=field.x;
    [element,exp]=sort(field);
    exp=double(exp);
    for i=1:16  
        [,col,]=find(H(i,:));
        shift=exp(nonzeros(H(i,:))');
        %Baux_2=17+(-(15-mod(shift,15)+1));
        Baux_2=(0-(15-mod(shift,15)+1));
        Qmn_prime((i-1)*16+1,col)=Qmn((i-1)*16+1,col);
        %Rmn_prime_aux((i-1)*16+1,col)=Rmn((i-1)*16+1,col);
        for b=1: length(col)
            if (Baux_2(b)==1)
               % Qmn_prime((i-1)*16+2:(i-1)*16+16,col(b))=[Qmn((i-1)*16+3:(i-1)*16+16,col(b)); Qmn((i-1)*16+2,col(b))];
               Qmn_prime((i-1)*16+2:(i-1)*16+16,col(b))=[Qmn((i-1)*16+3:(i-1)*16+16,col(b)); Qmn((i-1)*16+2,col(b))];
            else
                %Qmn_prime((i-1)*16+2:(i-1)*16+16,col(b))=[Qmn((i-1)*16+(Baux_2(b)):(i-1)*16+16,col(b)); Qmn((i-1)*16+2:(i-1)*16+Baux_2(b)-1,col(b)) ];
                Qmn_prime((i-1)*16+2:(i-1)*16+16,col(b))=[Qmn((i-1)*16+2-(Baux_2(b)):(i-1)*16+16,col(b)); Qmn((i-1)*16+2:(i-1)*16+2-Baux_2(b)-1,col(b)) ];
            end 
%             Rmn_prime_aux((i-1)*16+2:(i-1)*16+16,col(b))=circshift(Rmn((i-1)*16+2:(i-1)*16+16,col(b)),[Baux_2(b), 0]);
%             if(nonzeros(Rmn_prime_aux~=Rmn_prime))
%                 stop=1;
%             end
        end
        clear col
    end
end