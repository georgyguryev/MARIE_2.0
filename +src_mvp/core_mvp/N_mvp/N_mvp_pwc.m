function Jout = N_mvp_pwc(Jin, fN, dom_dims, op_dims)

% get dimensions of input current and N operator
L = dom_dims(1);
M = dom_dims(2);
N = dom_dims(3);

% reshape/allocate currents
Jin_tensor = reshape(Jin, dom_dims);
Jout       = Jin_tensor;

% apply N operator to x, y, z components
fJin_q = fftn(Jin_tensor(:,:,:,1),op_dims(1:3));

Jout1 = fN(:,:,:,1) .* fJin_q;
Jout2 = fN(:,:,:,2) .* fJin_q;
Jout3 = fN(:,:,:,3) .* fJin_q;

fJin_q = fftn(Jin_tensor(:,:,:,2),op_dims(1:3));
% y component of JIn, add contribution on 3 components of Jout
Jout1 = Jout1 + fN(:,:,:,2) .* fJin_q;
Jout2 = Jout2 + fN(:,:,:,4) .* fJin_q;
Jout3 = Jout3 + fN(:,:,:,5) .* fJin_q;

fJin_q = fftn(Jin_tensor(:,:,:,3),op_dims(1:3));
% z component of JIn, add contribution on 3 components of Jout
Jout1 = Jout1 + fN(:,:,:,3) .* fJin_q;
Jout2 = Jout2 + fN(:,:,:,5) .* fJin_q;
Jout3 = Jout3 + fN(:,:,:,6) .* fJin_q;

% apply ifft
Jout1 = ifftn(Jout1);
Jout2 = ifftn(Jout2);
Jout3 = ifftn(Jout3);

% cut noise, return relting fields
Jout(:,:,:,1) = Jout1(1:L,1:M,1:N);
Jout(:,:,:,2) = Jout2(1:L,1:M,1:N);
Jout(:,:,:,3) = Jout3(1:L,1:M,1:N);

% reshape Jout to vector; return result
Jout = reshape(Jout, size(Jin));