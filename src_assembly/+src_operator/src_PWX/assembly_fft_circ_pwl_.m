function [fA] = assembly_fft_circ_pwl_(A)

    [n1,n2,n3,n4,n5] = size(A);    

    if n5 == 6
        G_n1   = [ +1 -1 -1 +1 +1 +1 ];
        G_n2   = [ +1 -1 +1 +1 -1 +1 ];
        G_n3   = [ +1 +1 -1 +1 -1 +1 ];
        G_n1n2  = [ +1 +1 -1 +1 -1 +1 ];
        G_n1n3  = [ +1 -1 +1 +1 -1 +1 ];
        G_n2n3  = [ +1 -1 -1 +1 +1 +1 ];
        G_n1n2n3 = [ +1 +1 +1 +1 +1 +1 ];
    elseif n5 == 3
        G_n1   = [ -1 +1 +1 ];
        G_n2   = [ +1 -1 +1 ];
        G_n3   = [ +1 +1 -1 ];
        G_n1n2  = [ -1 -1 +1 ];
        G_n1n3  = [ -1 +1 -1 ];
        G_n2n3  = [ +1 -1 -1 ];
        G_n1n2n3 = [ -1 -1 -1 ];
    end

    TB_n1   = [ +1 -1 +1 +1 ];
    TB_n2   = [ +1 +1 -1 +1 ];
    TB_n3   = [ +1 +1 +1 -1 ];
    TB_n1n2  = [ +1 -1 -1 +1 ];
    TB_n1n3  = [ +1 -1 +1 -1 ];
    TB_n2n3  = [ +1 +1 -1 -1 ];
    TB_n1n2n3 = [ +1 -1 -1 -1 ];
    
    fA = zeros(2*n1,2*n2,2*n3,n4,n5);
    fA(1:n1,1:n2,1:n3,:,:) = A;

    ind_l  =  [1 2 3 4 1 1 1 2 2 3];
    ind_lp = [1 2 3 4 2 3 4 3 4 4];

    for j = 1:n4
        l = ind_l(j);
        lp = ind_lp(j);

        for i = 1:n5
            fA(n1+2:2*n1, 1:n2          , 1:n3, j, i)          = A(n1:-1:2 , 1:n2     , 1:n3, j, i)       * G_n1(i)        * TB_n1(l)        * TB_n1(lp);
            fA(1:n1        , n2+2:2*n2 , 1:n3, j, i)           = A(1:n1     , n2:-1:2 , 1:n3, j, i)       * G_n2(i)        * TB_n2(l)        * TB_n2(lp);
            fA(1:n1        , 1:n2          , n3+2:2*n3, j, i)  = A(1:n1     , 1:n2     , n3:-1:2, j, i)   * G_n3(i)        * TB_n3(l)        * TB_n3(lp);
            fA(n1+2:2*n1, n2+2:2*n2 , 1:n3, j, i)          = A(n1:-1:2 , n2:-1:2 , 1:n3, j, i)       * G_n1n2(i)     * TB_n1n2(l)    * TB_n1n2(lp);
            fA(n1+2:2*n1, 1:n2         , n3+2:2*n3, j, i)  = A(n1:-1:2 , 1:n2     , n3:-1:2, j, i)   * G_n1n3(i)     * TB_n1n3(l)    * TB_n1n3(lp);
            fA(1:n1        , n2+2:2*n2 , n3+2:2*n3, j, i)  = A(1:n1     , n2:-1:2 , n3:-1:2, j, i)   * G_n2n3(i)     * TB_n2n3(l)    * TB_n2n3(lp);
            fA(n1+2:2*n1, n2+2:2*n2 , n3+2:2*n3, j, i) = A(n1:-1:2  , n2:-1:2 , n3:-1:2, j, i)   * G_n1n2n3(i) * TB_n1n2n3(l) * TB_n1n2n3(lp);
        end
    end

    for j = 1:n4
            for i = 1:n5
                fA(:,:,:,j,i) = fftn(fA(:,:,:,j,i));
            end
    end

end
