function [Q]=L2Q_ver_2(reliability_mtrx2,H)
    Q=zeros(16*16,32);
    field=[gf(2,4).^(0:14)];
    field=field.x;
    [element,exp]=sort(field);
    exp=double(exp);
    for i=1:16  
        [,col,]=find(H(i,:));
        shift=exp(nonzeros(H(i,:))');
        Baux_2=(0-(15-mod(shift,15)+1));
        Q((i-1)*16+1,col)=reliability_mtrx2(1,col);
        for b=1: length(col)
            if (Baux_2(b)==1)
                Q((i-1)*16+2:(i-1)*16+16,col(b))=[reliability_mtrx2(16,col(b)); reliability_mtrx2(3-(Baux_2(b)):15,col(b))];
            else
                Q((i-1)*16+2:(i-1)*16+16,col(b))=[reliability_mtrx2(2-(Baux_2(b)):16,col(b)); reliability_mtrx2(2:2-Baux_2(b)-1,col(b)) ];
            end
        end 
        clear col
    end   