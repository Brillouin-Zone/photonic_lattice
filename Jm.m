% JD
fname = 'C:\von_Server\ETH\BSc Physics\7\Bachelorarbeit\plots';
k =1.5;
f = 4;
U = 1;
J_set = linspace(0,10,1000);
%D_set = linspace(0,10,numel(J_set));
m_set = linspace(0,10,numel(J_set)); % m = 1+2D / U => D = (m-1)*U /2


valid_sol_counter_Jm = zeros(numel(J_set), numel(J_set));

for Ji = 1: numel(J_set)
    J = J_set(Ji);
    for mi = 1:numel(m_set)
        m = m_set(mi);
        p_Jm = [1, -(2.*(((m-1)/2)+J)./U), ((((m-1)/2)+J).^2 + k^2/4)/U^2, -(f/U).^2];
        n0_Jm = roots(p_Jm);
        n_set1_Jm(Ji) = n0_Jm(1);
        n_set2_Jm(Ji) = n0_Jm(2);
        n_set3_Jm(Ji) = n0_Jm(3);
        if real(n_set1_Jm(Ji)) > 0 && imag(n_set1_Jm(Ji)) == 0
           valid_sol_counter_Jm(Ji, mi) = valid_sol_counter_Jm(Ji, mi) +1;
        end
        if real(n_set2_Jm(Ji)) > 0 && imag(n_set2_Jm(Ji)) == 0
           valid_sol_counter_Jm(Ji, mi) = valid_sol_counter_Jm(Ji, mi) +1;
        end
        if real(n_set3_Jm(Ji)) > 0 && imag(n_set3_Jm(Ji)) == 0
           valid_sol_counter_Jm(Ji, mi) = valid_sol_counter_Jm(Ji, mi) +1;
        end
    end
end
valid_sol_counter_Jm = flipud(rot90(valid_sol_counter_Jm));

mymap = [0 1 1
    1 1 0]; % cyan, yellow

figure
H_Jm = imagesc((valid_sol_counter_Jm));
colormap(mymap);
xt_Jm = get(gca, 'XTick');     
set(gca, 'XTick', xt_Jm, 'XTickLabel', xt_Jm/numel(J_set) * m_set(length(m_set))) ; 
yt_Jm = get(gca, 'YTick');     
set(gca, 'YTick', yt_Jm, 'YTickLabel', yt_Jm/numel(J_set)* J_set(length(J_set))) ; 
xlabel('m = 1+2\Delta / U');
ylabel('J/U');
%title('spinodals and critical point');
saveas(gcf, fullfile(fname, 'Jm.eps'), 'epsc'); 
saveas(gcf, 'Jm.pdf'); 
