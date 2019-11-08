% spinodals in (x, y) = (k, f)
% MFT: J = 1; U = const
U = 1;
k_set_kf = linspace(0, 2, 1000);
D = 1.5;
f_set_kf = linspace(0, 2, numel(k_set_kf));

fname = 'C:\von_Server\ETH\BSc Physics\7\Bachelorarbeit\plots';
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
ind_kf = find(sum(spinodal_2_kf,1)==0) ;
spinodal_2_kf(:,ind_kf) = [] ;


spinodal_1_k =  spinodal_1_kf(2, :) /numel(k_set_kf) * k_set_kf(length(k_set_kf));
spinodal_1_f =  spinodal_1_kf(1, :) /numel(f_set_kf) * f_set_kf(length(f_set_kf));
spinodal_2_k =  spinodal_2_kf(2, :) /numel(k_set_kf) * k_set_kf(length(k_set_kf));
spinodal_2_f =  spinodal_2_kf(1, :) /numel(f_set_kf) * f_set_kf(length(f_set_kf));
[xData1k, yData1f] = prepareCurveData( spinodal_1_k, spinodal_1_f );
[xData2k, yData2f] = prepareCurveData( spinodal_2_k, spinodal_2_f );
ft_kf = fittype( 'poly3' );
[fitresult1kf, gof1kf] = fit( xData1k, yData1f, ft_kf );
[fitresult2kf, gof2kf] = fit( xData2k, yData2f, ft_kf );

coeff_1kf = coeffvalues(fitresult1kf);
    %Spinodale oben: 0.0188    0.0530    0.0223    0.7044
coeff_2kf = coeffvalues(fitresult2kf);
    %Spinodale unten: -0.0182    0.0162    0.6050    0.0020
kappa_array = linspace(0, 2, 1000);
spinodale_1_fit = coeff_1kf(1)*kappa_array(:).^3 + coeff_1kf(2)*kappa_array(:).^2 + coeff_1kf(3)*kappa_array(:) + coeff_1kf(4);
spinodale_2_fit = coeff_2kf(1)*kappa_array(:).^3 + coeff_2kf(2)*kappa_array(:).^2 + coeff_2kf(3)*kappa_array(:) + coeff_2kf(4);

figure
title('fit of spinodals')
hold on
plot(kappa_array(:), spinodale_1_fit(:), 'r-', 'Linewidth', 0.8);
hold on
plot(kappa_array(:), spinodale_2_fit(:), 'r-', 'Linewidth', 0.8);
xlabel('\kappa / U');
ylabel('f / U');
%legend('f/U \approx 0.02(^{\kappa}/_{U})^3 + 0.05(^{\kappa}/_{U})^2 +0.02(^{\kappa}/_{U})+0.70(^{\kappa}/_{U})', 'f/U = -0.02(^{\kappa}/_{U})^3 + 0.02(^{\kappa}/_{U})^2 +0.61(^{\kappa}/_{U})+0.00(^{\kappa}/_{U})')
%XARROW1 = [0.5 0.56];
%YARROW1 = [0.8 0.7];
%XARROW2 = [0.45 0.35];
%YARROW2 = [0.25 0.35];
%annotation('textarrow',XARROW1,YARROW1,'String','f/U \approx 0.02(\kappa / U)^3 + 0.05(\kappa / U)^2 +0.02(\kappa / U)+0.70(\kappa / U) ')
%annotation('textarrow',XARROW2,YARROW2,'String','f/U \approx -0.02(\kappa / U)^3 + 0.02(\kappa / U)^2 +0.61(\kappa / U)+0.00(\kappa / U) ')


figure
H_kf = imagesc((valid_sol_counter_kf));
%colorbar
hold on
plot(spinodal_1_kf(2, length(spinodal_1_kf)), spinodal_1_kf(1, length(spinodal_1_kf)), 'k.', 'MarkerSize', 15) % critical point: last
%column of spinodal_1
hold on
%plot(kappa_array(:), spinodale_1_fit(:), 'r-');
plot(spinodal_1_kf(2, :), spinodal_1_kf(1, :), 'r-', 'Linewidth', 0.8);
hold on
plot(spinodal_2_kf(2, :), spinodal_2_kf(1, :), 'r-', 'Linewidth', 0.8);
%set(Plot, 'EdgeColor','none', 'FaceColor','interp');
%alpha(.1);
set(gca, 'YDIR', 'normal');
xt_kf = get(gca, 'XTick');     
set(gca, 'XTick', xt_kf, 'XTickLabel', xt_kf/numel(k_set_kf) * k_set_kf(length(k_set_kf))) ; 
yt_kf = get(gca, 'YTick'); 
set(gca, 'YTick', yt_kf, 'YTickLabel', yt_kf/numel(k_set_kf) * f_set_kf(length(f_set_kf)));  
ylabel('f / U');
xlabel('\kappa / U');
%title('spinodals and critical point');
%saveas(gcf, fullfile(fname, 'kf.eps'), 'epsc'); 
%saveas(gcf, 'kf.pdf'); 



figure
plot(spinodal_1_kf(2, length(spinodal_1_kf)), spinodal_1_kf(1, length(spinodal_1_kf)), 'k.', 'MarkerSize', 15) % critical point: last
%column of spinodal_1
hold on
plot(spinodal_1_kf(2, :), spinodal_1_kf(1, :), 'r-', 'Linewidth', 0.8);
%plot(kappa_array(:), spinodale_1_fit(:), 'r-');
hold on
plot(spinodal_2_kf(2, :), spinodal_2_kf(1, :), 'r-', 'Linewidth', 0.8);
%plot(kappa_array(:), spinodale_2_fit(:), 'g-');
set(gca, 'YDIR', 'normal');
xt_kf = get(gca, 'XTick');     
set(gca, 'XTick', xt_kf, 'XTickLabel', xt_kf/numel(k_set_kf) * k_set_kf(length(k_set_kf))) ; 
yt_kf = get(gca, 'YTick'); 
set(gca, 'YTick', yt_kf, 'YTickLabel', yt_kf/numel(k_set_kf) * f_set_kf(length(f_set_kf)));
dim1 = [.3 .6 .3 .3];
dim2 = [.3 .01 .3 .3];
dim3 = [.3 .3 .3 .3];
str1 = {'high-density solution','bright'};
str2 = {'low-density solution','dark'};
str3 = 'bistable';
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
annotation('textbox',dim2,'String',str2,'FitBoxToText','on');
annotation('textbox',dim3,'String',str3,'FitBoxToText','on');
xarrow = [0.5 0.5];
yarrow = [0.35 0.8];
annotation('textarrow',xarrow,yarrow);
hold on
plot(429, 387, '*g', 'MarkerSize', 12);
plot(429, 260, 'og', 'MarkerSize', 12);
ylabel('f / U');
xlabel('\kappa / U');
%title('spinodals and critical point');
saveas(gcf, fullfile(fname, 'kf_brightness.eps'), 'epsc'); 
saveas(gcf, 'kf_brightness.pdf'); 
%}
