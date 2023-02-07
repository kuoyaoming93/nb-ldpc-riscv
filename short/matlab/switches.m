function [wires]=switches()
    wires=zeros(2,16,16);
    for i=1:15
        for a=0:14
            for b=0:14
                if ((i-1)==a)
                    wires(:,a+2,i+1)=[i-1 -1]';
                end
                if((gf(2,4)^a+gf(2,4)^b)==gf(2,4)^(i-1))
                    wires(:,a+2,i+1)=[a b]';
                end
            end
        end
        wires(:,1,i+1)=[-1 i-1]';
    end
    wires(:,:,1)=[-1:14; -1:14];
end
