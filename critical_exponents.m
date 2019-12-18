%%% Simplification and derivation of the critical exponents in the rh-frame
 syms U nC h r delta_n
 
 Delta = 3*U*nC / 2;
 kappaC = sqrt(3) * U * nC;
 fC = U*nC^(3/2);
 n = nC * (delta_n +1);
    f_r0 = (4*h / sqrt(3) + sqrt(2/3)*0)/6 + fC;
    kappa_r0 = (4*0 / sqrt(3) - sqrt(2/3)*h)/6 + kappaC;
    f_h0 = (4*0 / sqrt(3) + sqrt(2/3)*r)/6 + fC;
    kappa_h0 = (4*r / sqrt(3) - sqrt(2/3)*0)/6 + kappaC;
 
 cubic_r0 = n^3 - 2*Delta*n^2 / U + (Delta^2 + kappa_r0^2 /4)*(n/U^2) - f_r0^2/U^2 == 0;
 cubic_h0 = n^3 - 2*Delta*n^2 / U + (Delta^2 + kappa_h0^2 /4)*(n/U^2) - f_h0^2/U^2 == 0;
 Ecubic_r0 = expand(simplify(cubic_r0));
 Ecubic_h0 = simplify(expand(cubic_h0));
    
     
 % for delta_n small, we neglegt the cubic term,
 Ecubic_r0_dn3 = subs(Ecubic_r0, delta_n^3, 0);  % equation between (41) and (42)
 Ecubic_h0_dn3 = subs(Ecubic_h0, delta_n^3, 0);  % equation between (42) and (43) 
     solution_r0_dn3 = solve(Ecubic_r0_dn3, delta_n);   % equation (42)
     solution_h0_dn3 = solve(Ecubic_h0_dn3, delta_n);   % equation (43)
            Taylor_solution_r0_dn3 = taylor(solution_r0_dn3, h, 'ExpansionPoint',0,'Order',2); % Taylor expansion at h=0 up to third order
                % RESULT:
                % delta n = - h*((18*2^(1/2)*U*nC^2 + 96*3^(1/2)*U*nC^(3/2))/(648*U^2*nC^3) - (2^(1/2)*(nC - 32))/(36*U*nC^2))
                %           - (2^(1/2)*(18*2^(1/2)*U*nC^2 + 96*3^(1/2)*U*nC^(3/2)))/(36*U*nC^2)
            Taylor_solution_h0_dn3 = taylor(solution_h0_dn3, r, 'ExpansionPoint',0,'Order',2); % Taylor expansion at r=0
                % RESULT:
                % delta n = - r*((2*nC - 1)/(18*U*nC^2) - (18*U*nC^2 - 6*6^(1/2)*U*nC^(3/2))/(162*U^2*nC^3))
                %           - (18*U*nC^2 - 6*6^(1/2)*U*nC^(3/2))/(18*U*nC^2)
 
 
 % Check: for r, h =0, we should find delta_n =0
 Ecubic_r0_dn3_h0 = subs(Ecubic_r0_dn3, h, 0);  
 Ecubic_h0_dn3_r0 = subs(Ecubic_h0_dn3, r, 0);
    solution_rh0 = solve(Ecubic_r0_dn3_h0, delta_n); % result: delta_n = 0
    solution_hr0 = solve(Ecubic_h0_dn3_r0, delta_n); % result: delta_n = 0
 
 
 
