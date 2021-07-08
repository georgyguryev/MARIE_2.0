function Jout = K_mvp_pwc(Jin, fK, dom_dims, op_dims)


% get dimensions of input current and N operator
L = dom_dims(1);
M = dom_dims(2);
N = dom_dims(3);

% reshape/allocate currents
Jin_tensor = reshape(Jin, dom_dims);
Jout       = zeros(dom_dims);


% ---------------------------------------------------------------------
% apply fft and mv-op for each of the components of JIn
% ---------------------------------------------------------------------
   
% apply fft and mv-op
fJin_q = fftn(Jin_tensor(:,:,:,1), op_dims(1:3));
Jout3 =  fK(:,:,:,2) .* fJin_q;
Jout2 = -fK(:,:,:,3) .* fJin_q;

fJin_q = fftn(Jin_tensor(:,:,:,2), op_dims(1:3));
Jout1 = fK(:,:,:,3) .* fJin_q;
Jout3 = Jout3 - fK(:,:,:,1) .* fJin_q;

fJin_q = fftn(Jin_tensor(:,:,:,3), op_dims(1:3));
Jout2 = Jout2 + fK(:,:,:,1) .* fJin_q;
Jout1 = Jout1 - fK(:,:,:,2) .* fJin_q;

% apply ifft
Jout1 = ifftn(Jout1);
Jout2 = ifftn(Jout2);
Jout3 = ifftn(Jout3);

% apply 
Jout(:,:,:,1) = Jout1(1:L,1:M,1:N);
Jout(:,:,:,2) = Jout2(1:L,1:M,1:N);
Jout(:,:,:,3) = Jout3(1:L,1:M,1:N);

% reshape Jout to vector; return result
Jout = - reshape(Jout, size(Jin));