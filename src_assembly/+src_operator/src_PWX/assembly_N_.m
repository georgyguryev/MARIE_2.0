function I_V = assembly_N_(r,EM_const,res,pwx,L)

    [n1,n2,n3,~] = size(r);
    k0 = EM_const.k0;
    r_n = [0 0 0]';
    d = [res,res,res];

    Np_1D_far_V = 2;
    Np_1D_medium_V = 5;
    Np_1D_near_S = 12;
    
    I_V = zeros(n1,n2,n3,pwx,6);
    

    %% far distance cells - calculate the volume-volume integrals

    [W6D,X,Y,Z,Xp,Yp,Zp] = weights_points(Np_1D_far_V,6);

    parfor mx = 1:n1
        for my = 1:n2
            for mz = 1:n3

                m = [mx,my,mz];
                r_m = ((m-1) .* d)';
                I_V(mx,my,mz,:,:) = VV_Nop(W6D,X,Y,Z,Xp,Yp,Zp,k0,r_m,r_n,res,pwx);

            end
        end
    end

    %% medium distance cells - calculate the volume-volume integrals

    [W6D,X,Y,Z,Xp,Yp,Zp] = weights_points(Np_1D_medium_V,6);

    parfor mx = 1:Np_1D_medium_V
        for my = 1:Np_1D_medium_V
            for mz = 1:Np_1D_medium_V

                m = [mx,my,mz];
                r_m = ((m-1) .* d)';
                I_V(mx,my,mz,:,:) = VV_Nop(W6D,X,Y,Z,Xp,Yp,Zp,k0,r_m,r_n,res,pwx);

            end
        end
    end

    %% nearby cells - calculate the surface-surface integrals

    [W4D,A,B,C,D] = weights_points(Np_1D_near_S,4);
    [I1_co,I2_co,I3_co,I4_co] = surface_surface_coeff_Nop(res,L);

    Np_s_1 = 2;
    Np_s_2 = 2;
    Np_s_3 = 2;
    if n1 <2
        Np_s_1 = 1;
    end
    if n2 <2
        Np_s_2 = 1;
    end
    if n3 <2
        Np_s_3 = 1;
    end

    parfor mx = 1:Np_s_1
        for my = 1:Np_s_2
            for mz = 1:Np_s_3

                m = [mx,my,mz];
                r_m = ((m-1) .* d)';      
                I_V(mx,my,mz,:,:)  = surface_surface_kernels_Nop(I1_co,I2_co,I3_co,I4_co,W4D,A,B,C,D,Np_1D_near_S,k0,r_m,r_n,m,res,L);
           end
        end
    end

end
