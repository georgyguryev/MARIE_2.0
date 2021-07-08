function Jout = N_mvp_pwl_tucker_sparsefft(Jin,T,dom_dims,op_dims)

    n1 = dom_dims(1);
    n2 = dom_dims(2);
    n3 = dom_dims(3);
    pwx = dom_dims(end);
    m1 = op_dims(1);
    m2 = op_dims(2);
    m3 = op_dims(3);
    
    Jin_tensor = reshape(Jin, dom_dims);
    
    fJ  = zeros(m1,m2,n3,pwx,'like',Jin_tensor);
    NfJ = zeros(m1,m2,m3,'like',Jin_tensor);
    Jout  = zeros(n1,n2,n3,pwx,'like',Jin_tensor);

    pq = [1 2 3;2 4 5;3 5 6];
    llp = [1 5 6 7;-5 2 8 9;-6 8 3 10;-7 9 10 4];
    
    for i=1:pwx
        fJ(:,:,:,i) = fft(fft(Jin_tensor(:,:,:,i),m1,1),m2,2);
    end
    
    for p = 1:3
        for l = 1:4
            
            left = 4*(p-1) + l;
            
            for q = 1:3
                for lp = 1:4
                        
                    i = abs(llp(l,lp));
                    j = pq(p,q);
                    
                    right = 4*(q-1) + lp;
                    
                    N_Tucker = nmp(nmp(nmp(T(i,j).G,T(i,j).U1,1),T(i,j).U2,2),T(i,j).U3,3);

                    a = fft(fJ(1:n1,1:n2,:,right),m3,3);
                    NfJ(1:n1,1:n2,:) = NfJ(1:n1,1:n2,:) + sign(llp(l,lp)) * N_Tucker(1:n1,1:n2,:).*a;

                    a = fft(fJ(1:n1,n2+1:m2,:,right),m3,3);
                    NfJ(1:n1,n2+1:m2,:) = NfJ(1:n1,n2+1:m2,:) + sign(llp(l,lp)) * N_Tucker(1:n1,n2+1:m2,:).*a;

                    a = fft(fJ(n1+1:m1,1:n2,:,right),m3,3);
                    NfJ(n1+1:m1,1:n2,:) = NfJ(n1+1:m1,1:n2,:) + sign(llp(l,lp)) * N_Tucker(n1+1:m1,1:n2,:).*a;
                    
                    a = fft(fJ(n1+1:m1,n2+1:m2,:,right),m3,3);
                    NfJ(n1+1:m1,n2+1:m2,:)  =  NfJ(n1+1:m1,n2+1:m2,:) + sign(llp(l,lp)) * N_Tucker(n1+1:m1,n2+1:m2,:).*a;
                    
                end
            end
            
            clear a N_Tucker;
            
            NfJ = ifftn(NfJ);
            Jout(:,:,:,left) = NfJ(1:n1,1:n2,1:n3);
            NfJ = zeros(m1,m2,m3,'like',Jin);
            
        end
    end
    
    Jout = reshape(Jout, size(Jin));

end