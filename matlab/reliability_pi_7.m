function [L_n_aux,L_n]=reliability_pi_7(rx_points,sigma)
    q_field=16;
    m_field=4;
    N_code=32;
    comp=gf(zeros(1,(q_field-1)));
    comp(2:q_field)=gf(2,m_field).^(0:(q_field-2));
    comp=dec2bin(comp.x)*1-48;
    comp=1-2*comp;
    for i=1:N_code  
        for i2 =1:q_field
            d=1; 
            d=(1/sigma^2)*rx_points(i,:)/2;
            L_n_aux(i2,i)=sum(d.*comp(i2,:)); 
        end  
         L_n(:,i)=abs(L_n_aux(:,i)-max(L_n_aux(:,i)));
    end 