function Jout = K_mvp_pwl(Jin, T, dom_dims, op_dims)

n1 = dom_dims(1);
n2 = dom_dims(2);
n3 = dom_dims(3);
pwx = dom_dims(end);
m1 = op_dims(1);
m2 = op_dims(2);
m3 = op_dims(3);

Jin = reshape(Jin, n1,n2,n3,pwx);
fJ  = zeros(m1,m2,m3,pwx,'like',Jin);
KfJ = zeros(m1,m2,m3,'like',Jin);
Jout  = zeros(n1,n2,n3,pwx,'like',Jin);

pq = [0 -3 +2;+3 0 -1;-2 +1 0];
llp = [1 5 6 7;-5 2 8 9;-6 8 3 10;-7 9 10 4];

for i=1:pwx
    fJ(:,:,:,i) = fftn(Jin(:,:,:,i),[m1,m2,m3]);
end

q_v = 1:3;

for p = 1:3
    for l = 1:4
        
        left = 4*(p-1) + l;
        
        for q = q_v(1:end ~=p)
            for lp = 1:4
                
                right = 4*(q-1) + lp;
                KfJ = KfJ + sign(pq(p,q)) * sign(llp(l,lp)) * T(:,:,:,abs(llp(l,lp)),abs(pq(p,q))) .* fJ(:,:,:,right);
                
            end
        end
        
        KfJ = ifftn(KfJ);
        Jout(:,:,:,left) = KfJ(1:n1,1:n2,1:n3);
        KfJ = zeros(m1,m2,m3,'like',Jin);
        
    end
end


% reshape Jout to vector; return result
Jout = reshape(Jout, size(Jin));