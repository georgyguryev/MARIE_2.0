function Jout = K_mvp_pwc_tucker(Jin,T,dom_dims,op_dims)
    
n1 = dom_dims(1);
n2 = dom_dims(2);
n3 = dom_dims(3);
pwx = dom_dims(end);
m1 = op_dims(1);
m2 = op_dims(2);
m3 = op_dims(3);

Jin = reshape(Jin,n1,n2,n3,pwx);
fJ  = zeros(m1,m2,m3,pwx,'like',Jin);
KfJ = zeros(m1,m2,m3,'like',Jin);
Jout  = zeros(n1,n2,n3,pwx,'like',Jin);

pq = [0 -3 +2;+3 0 -1;-2 +1 0];

for i=1:pwx
    fJ(:,:,:,i) = fftn(Jin(:,:,:,i),[m1,m2,m3]);
end

q_v = 1:3;

for p = 1:3
    
    for q = q_v(1:end ~=p)
        
        j = abs(pq(p,q));
        
        KfJ = KfJ + sign(pq(p,q)) * nmp(nmp(nmp(T(j).G,T(j).U1,1),T(j).U2,2),T(j).U3,3) .* squeeze(fJ(:,:,:,q));
        
    end
    
    KfJ = ifftn(KfJ);
    Jout(:,:,:,p) = KfJ(1:n1,1:n2,1:n3);
    KfJ = zeros(m1,m2,m3,'like',Jin);
    
end




% reshape Jout to vector; return result
Jout = reshape(Jout, size(Jin));