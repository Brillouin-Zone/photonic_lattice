% spinodals in (x, y) = (k, f)
% MFT: J = 1; U = const
U = 1;
k_set_kf = linspace(0, 2, 1000);
D = 1.5;
f_set_kf = linspace(0, 5, numel(k_set_kf));

fname = 'C:\von_Server\ETH\BSc Physics\7\Bachelorarbeit\MATLAB';
counter_kf = 1;
valid_sol_counter_kf = zeros(numel(k_set_kf), numel(f_set_kf)); 
spinodal_1_kf = zeros(2, numel(k_set_kf));
spinodal_2_kf = zeros(2, numel(k_set_kf));
for fi = 1: numel(f_set_kf)
    f = f_set_kf(fi);
    for ki = 1:numel(k_set_kf)
        k = k_set_kf(ki);
        p_kf = [1, -(2.*D./U), ((D.^2 + k.^2./4)/U.^2), -(f./U)^2];
        n0_kf = roots(p_kf);
        n_set1_kf(fi) = n0_kf(1); % all complex
        n_set2_kf(fi) = n0_kf(2); % all complex
        n_set3_kf(fi) = n0_kf(3); % all real
        if real(n_set1_kf(fi)) > 0 && (imag(n_set1_kf(fi)) == 0 | abs(imag(n_set1_kf(fi))) <= 10^(-10) )
            valid_sol_counter_kf(ki, fi) = valid_sol_counter_kf(ki, fi) +1;
        end
        if real(n_set2_kf(fi)) > 0 && (imag(n_set2_kf(fi)) == 0 | abs(imag(n_set2_kf(fi))) <= 10^(-10) )
                valid_sol_counter_kf(ki, fi) = valid_sol_counter_kf(ki, fi) +1;
        end
        if real(n_set3_kf(fi)) > 0 && (imag(n_set3_kf(fi)) == 0 | abs(imag(n_set3_kf(fi))) <= 10^(-10) )
                    valid_sol_counter_kf(ki, fi) = valid_sol_counter_kf(ki, fi) +1;
        end
    end
end
valid_sol_counter_kf = flipud(rot90(valid_sol_counter_kf));
for i = 1: size(valid_sol_counter_kf , 1)
    for j = 1: size(valid_sol_counter_kf , 2)-1
        if valid_sol_counter_kf(i,j) ==1 && valid_sol_counter_kf(i, j+1) == 3 
            spinodal_1_kf(1, i) =  i;% k-coordinate of the corresponding points
            spinodal_1_kf(2, i) =  j;% f-coordinate of the corresponding points
        elseif j < size(spinodal_1_kf, 2) && spinodal_1_kf(1, j) == 0
            spinodal_1_kf(:, j) = []; % the coordinates of the critical point are in the last column
        end
        if valid_sol_counter_kf(i,j) == 3 && valid_sol_counter_kf(i, j+1) == 1
            spinodal_2_kf(1, i) =  i;% k-coordinate of the corresponding points
            spinodal_2_kf(2, i) =  j;% f-coordinate of the corresponding points
        elseif i < size(spinodal_2_kf, 2) && spinodal_2_kf(1, i) == 0
                spinodal_2_kf(:, i) = []; 
        end
    end
end
figure
H_kf = imagesc((valid_sol_counter_kf));
colorbar
hold on
plot(spinodal_1_kf(2, length(spinodal_1_kf)), spinodal_1_kf(1, length(spinodal_1_kf)), 'm*') % critical point: last
%column of spinodal_1
hold on
plot(spinodal_1_kf(2, :), spinodal_1_kf(1, :), 'r-')
hold on
plot(spinodal_2_kf(2, :), spinodal_2_kf(1, :), 'g-')
set(gca, 'YDIR', 'normal');
xt_kf = get(gca, 'XTick');     
set(gca, 'XTick', xt_kf, 'XTickLabel', xt_kf/numel(k_set_kf) * k_set_kf(length(k_set_kf))) ; 
yt_kf = get(gca, 'YTick'); 
set(gca, 'YTick', yt_kf, 'YTickLabel', yt_kf/numel(k_set_kf) * f_set_kf(length(f_set_kf)));  
ylabel('f / U');
xlabel('\kappa / U');
title('spinodals and critical point');
saveas(gcf, fullfile(fname, 'kf.eps'), 'epsc'); 
saveas(gcf, 'kf.pdf'); 

spinodal_1_k =  spinodal_1_kf(2, :) /numel(k_set_kf) * k_set_kf(length(k_set_kf));
spinodal_1_f =  spinodal_1_kf(1, :) /numel(f_set_kf) * f_set_kf(length(f_set_kf));
spinodal_2_k =  spinodal_2_kf(2, :) /numel(k_set_kf) * k_set_kf(length(k_set_kf));
spinodal_2_f =  spinodal_2_kf(1, :) /numel(f_set_kf) * f_set_kf(length(f_set_kf));
[xData1k, yData1f] = prepareCurveData( spinodal_1_k, spinodal_1_f );
[xData2k, yData2f] = prepareCurveData( spinodal_2_k, spinodal_2_f );
ft_kf = fittype( 'smoothingspline' );
[fitresult1kf, gof1kf] = fit( xData1k, yData1f, ft_kf );
[fitresult2kf, gof2kf] = fit( xData2k, yData2f, ft_kf );

