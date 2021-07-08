function H_tot = Jcb2Htot(Jc, Jb, K_ext, Zc2b_Kop, P,S)


Jtot      = P * Jc + S * Jb;
H_vox_ext = K_ext(Jtot);


H_tot = gather(Zc2b_Kop * Jc +  S.' * H_vox_ext(:));