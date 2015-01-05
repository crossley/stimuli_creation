function A = make_covar_matrix(d1,d2)

A = [d1+d2 d1-d2;d1-d2 d1+d2]/2;