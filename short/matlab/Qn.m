function [Q_n,index,mag]=Qn(H,Rmn_prime,reliability_mtrx2)
    Q_n=zeros(16,32);
    index=zeros(1,32);
    mag=zeros(1,32);
    for i=1:32
        [row,  ]=find(H(:,i));
        row=(row-1)*16+1;
        accum=zeros(16,1);
        for b=1:length(row)
            accum=Rmn_prime(row(b):(row(b)+15),i)+accum;
        end
        Q_n(:,i)=accum+reliability_mtrx2(:,i);
        [ index(i),mag(i)]=min(Q_n(:,i));
    end
end