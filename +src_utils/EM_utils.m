classdef EM_utils < handle
    
    % Borrowed from Jose's GMT project
	properties (Constant)
		mu0  = 4e-7 * pi;
		c0   = 299792458;
		e0   = 1/ src_utils.EM_utils.c0^2 / src_utils.EM_utils.mu0;
		eta0 = 3.767303134617706e+002;
	end
	properties (SetAccess = immutable)
		freq
		omega
		lambda0
		k0
		alpha
	end
	properties (Dependent)
		ce
		cm
	end
	methods
		function emc = EM_utils(freq)
			emc.freq    = freq;
			emc.omega   = 2*pi*freq;
			emc.lambda0 = emc.c0/freq;
			emc.k0      = 2*pi/emc.lambda0;
			emc.alpha   = 2*pi*42.58; % 7T
        end
        
		function ce = get.ce(emc)
			ce = 1j*emc.omega * src_utils.EM_utils.e0;
		end
		function cm = get.cm(emc)
			cm = 1j*emc.omega * src_utils.EM_utils.mu0;
		end
	end
end