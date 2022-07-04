function [ G,swap ] = NB_LDPC_generator_matrix( A )
% Eliminación de Gauss-Jordan

H = A;
[M N] = size(A);
swap = 1:N;

for m = 0:M-1
    for n = m:N-1
        for k = m:M-1
            if (H(M-k,N-n)~=0) 
                break;
            end
        end
        if (H(M-k,N-n)~=0)
            break
        end
    end
    if (H(M-k,N-n)==0)
        m = m-1;
        break;
    end 
    if (n~=m)
        tmp = H(:,N-m);   
        H(:,N-m) = H(:,N-n); 
        H(:,N-n) = tmp;
        tmp = swap(N-m);
        swap(N-m) = swap(N-n);
        swap(N-n) = tmp;
    end
    if (k~=m)
        H(M-m,:) = H(M-m,:)+H(M-k,:);
        H(M-m,:) = H(M-m,:)/H(M-m,N-m);
    end
    if (k==m && H(M-m,N-m)~=1)
        H(M-m,:) = H(M-m,:)/H(M-m,N-m);
    end
    for k = m+1:M-1
        if (H(M-k,N-m)~=0)
            H(M-k,:) = H(M-k,:)/H(M-k,N-m)+H(M-m,:)/H(M-m,N-m); 
        end
    end
end
for mi = 0:m-1 
    for k = mi+1:m
        if(k==mi+1)
            H(M-m+k-1,:)= H(M-m+k-1,:)/H(M-m+k-1,N-m+mi);
        end
        if (H(M-m+k,N-m+mi)~=0) 
            H(M-m+k,:) = H(M-m+k,:)/H(M-m+k,N-m+mi)+H(M-m+mi,:)/H(M-m+mi,N-m+mi);
        end
    end
end
H(M,:)=H(M,:)/H(M,N);
G = [eye(N-m-1); H(M-m:M,1:N-m-1)]';    

