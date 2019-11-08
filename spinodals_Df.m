% spinodals in (x, y) = (D, f)
% MFT: J = 1; U = const
U_Df = 1;
k_Df = 1;
D_set_Df = linspace(0, 10, 1000);
f_set_Df = linspace(0, 10, numel(D_set_Df));
fname = 'C:\von_Server\ETH\BSc Physics\7\Bachelorarbeit\plots';
counter_Df = 1;
valid_sol_counter_Df = zeros(numel(D_set_Df), numel(f_set_Df)); 
spinodal_1_Df = zeros(2, numel(D_set_Df));
spinodal_2_Df = zeros(2, numel(D_set_Df));
for fi = 1: numel(f_set_Df)
    f_Df = f_set_Df(fi);
    for Di = 1:numel(D_set_Df)
        D_Df = D_set_Df(Di);
        p_Df = [1, -(2.*D_Df./U_Df), ((D_Df.^2 + k_Df.^2./4)/U_Df.^2), -(f_Df./U_Df).^2];
        n0_Df = roots(p_Df);
        n_set1_Df(fi) = n0_Df(1);
        n_set2_Df(fi) = n0_Df(2);
        n_set3_Df(fi) = n0_Df(3);
        if real(n_set1_Df(fi)) > 0 && imag(n_set1_Df(fi)) == 0
           valid_sol_counter_Df(Di, fi) = valid_sol_counter_Df(Di, fi) +1;
        end
        if real(n_set2_Df(fi)) > 0 && imag(n_set2_Df(fi)) == 0
           valid_sol_counter_Df(Di, fi) = valid_sol_counter_Df(Di, fi) +1;
        end
        if real(n_set3_Df(fi)) > 0 && imag(n_set3_Df(fi)) == 0
           valid_sol_counter_Df(Di, fi) = valid_sol_counter_Df(Di, fi) +1;
        end
    end
end
valid_sol_counter_Df = flipud(rot90(valid_sol_counter_Df));
for i = 1: size(valid_sol_counter_Df , 1)-1
    for j = 1: size(valid_sol_counter_Df , 2)
        if valid_sol_counter_Df(i,j) ==1 && valid_sol_counter_Df(i+1, j) == 3 
            spinodal_1_Df(1, i) =  i;% k-coordinate of the corresponding points
            spinodal_1_Df(2, i) =  j;% f-coordinate of the corresponding points
        elseif j < size(spinodal_1_Df, 2) && spinodal_1_Df(1, j) == 0
            spinodal_1_Df(:, j) = []; % the coordinates of the critical point are in the last column
        end
        if valid_sol_counter_Df(i,j) == 3 && valid_sol_counter_Df(i+1, j) == 1
            spinodal_2_Df(1, i) =  i;% k-coordinate of the corresponding points
            spinodal_2_Df(2, i) =  j;% f-coordinate of the corresponding points
        end
    end
end
ind_Df = find(sum(spinodal_2_Df,1)==0) ;
spinodal_2_Df(:,ind_Df) = [] ;

figure
H_Df = imagesc((valid_sol_counter_Df));
%colorbar
hold on
plot(spinodal_1_Df(2, 1), spinodal_1_Df(1, 1), 'k.', 'MarkerSize', 15);
hold on
plot(spinodal_1_Df(2, :), spinodal_1_Df(1, :), 'r-', 'Linewidth', 0.8);
hold on
plot(spinodal_2_Df(2, :), spinodal_2_Df(1, :), 'r-', 'Linewidth', 0.8);
%set(H_Df, 'EdgeColor','none', 'FaceColor','interp');
%alpha(.1);
set(gca, 'YDIR', 'normal');
xt_Df = get(gca, 'XTick');     
set(gca, 'XTick', xt_Df, 'XTickLabel', xt_Df/numel(D_set_Df) * D_set_Df(length(D_set_Df))) ; 
yt_Df = get(gca, 'YTick');     
set(gca, 'YTick', yt_Df, 'YTickLabel', yt_Df/numel(D_set_Df)* f_set_Df(length(f_set_Df))) ; 
ylabel('f / U');
xlabel('\Delta / U');
%title('spinodals and critical point');
%saveas(gcf, fullfile(fname, 'Df.eps'), 'epsc'); 
%saveas(gcf, 'Df.pdf'); 


figure
plot(spinodal_1_Df(2, 1), spinodal_1_Df(1, 1), 'k.', 'MarkerSize', 12);
hold on
plot(spinodal_1_Df(2, :), spinodal_1_Df(1, :), 'r-', 'Linewidth', 0.8);
hold on
plot(spinodal_2_Df(2, :), spinodal_2_Df(1, :), 'r-', 'Linewidth', 0.8);
set(gca, 'YDIR', 'normal');
xt_Df = get(gca, 'XTick');     
set(gca, 'XTick', xt_Df, 'XTickLabel', xt_Df/numel(D_set_Df) * D_set_Df(length(D_set_Df))) ; 
yt_Df = get(gca, 'YTick');     
set(gca, 'YTick', yt_Df, 'YTickLabel', yt_Df/numel(D_set_Df)* f_set_Df(length(f_set_Df))) ; 
%dim1 = [.5 .6 .3 .3];
%dim2 = [.5 .0 .3 .3];
dim3 = [.7 .3 .3 .3];
%str1 = {'high-density solution','bright'};
%str2 = {'low-density solution; dark'};
str3 = 'bistable';
%annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
%annotation('textbox',dim2,'String',str2,'FitBoxToText','on');
annotation('textbox',dim3,'String',str3,'FitBoxToText','on');
xarrow = [0.55 0.55];
yarrow = [0.15 0.6];
annotation('textarrow',xarrow,yarrow);
hold on
plot(540, 484, '*g', 'MarkerSize', 12);
hold on
plot(540, 116, 'og', 'MarkerSize', 12);
ylabel('f / U');
xlabel('\Delta / U');
saveas(gcf, fullfile(fname, 'Df_brightness.eps'), 'epsc'); 
saveas(gcf, 'Df_brightness.pdf'); 


spinodal_1_D =  spinodal_1_Df(2, :) /numel(D_set_Df) * D_set_Df(length(D_set_Df));
spinodal_1_F =  spinodal_1_Df(1, :) /numel(f_set_Df) * f_set_Df(length(f_set_Df));
spinodal_2_D =  spinodal_2_Df(2, :) /numel(D_set_Df) * D_set_Df(length(D_set_Df));
spinodal_2_F =  spinodal_2_Df(1, :) /numel(f_set_Df) * f_set_Df(length(f_set_Df));
[xData1D, yData1F] = prepareCurveData( spinodal_1_D, spinodal_1_F );
[xData2D, yData2F] = prepareCurveData( spinodal_2_D, spinodal_2_F );
ft_Df = fittype( 'poly3' );
[fitresult1Df, gof1Df] = fit( xData1D, yData1F, ft_Df );
[fitresult2Df, gof2Df] = fit( xData2D, yData2F, ft_Df );

coeff_1Df = coeffvalues(fitresult1Df);
    %Spinodale oben: 0.0009   -0.0213    0.2570    0.2566
coeff_2Df = coeffvalues(fitresult2Df);
    %Spinodale unten: -0.0029    0.1124    0.3840   -0.0229
Delta_array = linspace(0, 10, 1000);
spinodale_1_Fit = coeff_1Df(1)*Delta_array(:).^3 + coeff_1Df(2)*Delta_array(:).^2 + coeff_1Df(3)*Delta_array(:) + coeff_1Df(4);
spinodale_2_Fit = coeff_2Df(1)*Delta_array(:).^3 + coeff_2Df(2)*Delta_array(:).^2 + coeff_2Df(3)*Delta_array(:) + coeff_2Df(4);

    figure
    title('fit of spinodals')
    plot(Delta_array(:), spinodale_1_Fit(:), 'r-');
    hold on
    plot(Delta_array(:), spinodale_2_Fit(:), 'g-');
    xlabel('\Delta / U');
    ylabel('f / U');








