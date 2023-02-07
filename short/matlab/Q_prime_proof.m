function [Q_prime]=Q_prime_proof(H,Rmn_prime,reliability_mtrx2)
    Q_prime=zeros(16*16,32);
    for i=1:32
        [row,  ]=find(H(:,i));
        row=(row-1)*16+1;
        for a=1:length(row)
            accum=zeros(16,1);
            for b=1:length(row)
                if a==b
                else
                    accum=Rmn_prime(row(b):(row(b)+15),i)+accum;
                end
            end
            Q_prime(row(a):(row(a)+15),i)=accum+reliability_mtrx2(:,i);
            Q_prime(row(a):(row(a)+15),i)=Q_prime(row(a):(row(a)+15),i)-min(Q_prime(row(a):(row(a)+15),i)); 
        end     
    end
end