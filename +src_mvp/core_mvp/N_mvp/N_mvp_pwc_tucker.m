function Jout = N_mvp_pwc_tucker(Jin,T,dom_dims,op_dims)

    n1 = dom_dims(1);
    n2 = dom_dims(2);
    n3 = dom_dims(3);
    pwx = dom_dims(end);
    m1 = op_dims(1);
    m2 = op_dims(2);
    m3 = op_dims(3);
    
    Jin_tensor = reshape(Jin, dom_dims);
    
    fJ  = zeros(m1,m2,m3,pwx,'like',Jin_tensor);
    NfJ = zeros(m1,m2,m3,'like',Jin_tensor);
    Jout  = zeros(n1,n2,n3,pwx,'like',Jin_tensor);

    pq = [1 2 3;2 4 5;3 5 6];

    for i=1:pwx
        fJ(:,:,:,i) = fftn(Jin_tensor(:,:,:,i),[m1,m2,m3]);
    end

    for p = 1:3

        for q = 1:3

            j = pq(p,q);
            NfJ  =  NfJ + nmp(nmp(nmp(T(j).G,T(j).U1,1),T(j).U2,2),T(j).U3,3) .* squeeze(fJ(:,:,:,q));

        end

        NfJ = ifftn(NfJ);
        Jout(:,:,:,p) = NfJ(1:n1,1:n2,1:n3);
        NfJ = zeros(m1,m2,m3,'like',Jin_tensor);

    end
    
    Jout = reshape(Jout, size(Jin));
        
end
