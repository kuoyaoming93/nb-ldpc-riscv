function [R_Forward,R_Backward,R_Aux,R,col]=R_proof(row,Q_2,H,wires)
    [,col,]=find(H(row,:)); 
    R=Q_2(((row-1)*16+1):(((row-1)*16)+16), col);
    %%Forward metrics
    %%Backward metrics
    R_Forward=zeros(16,4);
    R_Forward(:,1)=R(:,1);
    R_Backward=zeros(16,4);
    R_Backward(:,4)=R(:,4);
    for a=2:3
            R_aux=R(:,a);
            R_aux_2=R(:,5-a);
            R_Forward(:,a)=min(max(repmat(R_Forward(:,a-1),1,16),R_aux(wires+2)))';
            R_Backward(:,5-a)=min(max(repmat(R_Backward(:,5-a+1),1,16),R_aux_2(wires+2)))';
    end
    %%Standard Min-Max
    R_Aux=zeros(16,4);
    R_Aux(:,4)=R(:,4);
    for a=2:3
            R_Backward_aux=(R_Backward(:,a+1));
            R_Aux(:,a)=min(max(repmat(R_Forward(:,a-1),1,16),R_Backward_aux(wires+2)))';
    end
    R_Aux(:,1)=R_Backward(:,2);
    R_Aux(:,4)=R_Forward(:,3);
end   