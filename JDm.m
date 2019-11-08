% JDm
     % (x,y,z) = (J, D, m)
fname = 'C:\von_Server\ETH\BSc Physics\7\Bachelorarbeit\plots';
U = 3;
f = 1;
k = 1;
J_set = linspace(0, 4, 100);
D_set = linspace(0,4 ,numel(J_set));
DJ_set = (D_set(:) + J_set(:)).';
m_set = 1 +2*DJ_set(:) / U;
Valid_sol_counter_JDm = zeros(numel(J_set), numel(J_set), numel(J_set)); 
%critical_m = zeros(1, numel(J_set));
%critical_J = zeros(1, numel(J_set));
%critical_D = zeros(1, numel(J_set));

for Ji = 1:numel(J_set)
    J = J_set(Ji);
    for Di = 1:numel(J_set)
        D = D_set(Di);
        for mi = 1:numel(J_set)
            m = m_set(mi);
            P_JDm = [1, -(2.*(D+J)./U), (((D+J).^2 + k.^2./4)/U.^2), -(f./U)^2];
            N0_JDm = roots(P_JDm);
            N_set1_JDm(Ji, Di, mi) = N0_JDm(1);
            N_set2_JDm(Ji, Di, mi) = N0_JDm(2);
            N_set3_JDm(Ji, Di, mi) = N0_JDm(3);
            if real(N_set1_JDm(Ji, Di, mi)) > 0 && (imag(N_set1_JDm(Ji, Di, mi)) == 0)
                Valid_sol_counter_JDm(Ji, Di, mi) = Valid_sol_counter_JDm(Ji, Di, mi) +1;
            end
            if real(N_set2_JDm(Ji, Di, mi)) > 0 && (imag(N_set2_JDm(Ji, Di, mi)) == 0)
                Valid_sol_counter_JDm(Ji, Di, mi) = Valid_sol_counter_JDm(Ji, Di, mi) +1;
            end
            if real(N_set3_JDm(Ji, Di, mi)) > 0 && (imag(N_set3_JDm(Ji, Di, mi)) == 0)
                Valid_sol_counter_JDm(Ji, Di, mi) = Valid_sol_counter_JDm(Ji, Di, mi) +1;
            end
        end
    end
end
%Valid_sol_counter_JDm = rot90(fliplr(Valid_sol_counter_JDm));

%{
b = 0;
for Di = 1:numel(J_set) % kappa ebenen
    for Ji = 1:numel(J_set) % 
        for mi = 1:numel(J_set)
            if Valid_sol_counter_JDm(Di, Ji, mi) == 3
                critical_J(Ji) = Ji;
                critical_m(Ji) = mi;
                critical_D(Ji) = Di;
                b = 1;
                break;
            end
        end
        if b == 1
            b = 0;
            break;
        end
    end
end
%}

%Valid_sol_counter_JDm = rot90(fliplr(Valid_sol_counter_JDm));

mymap = [0 1 1
    1 1 0]; % cyan, yellow

figure
[mesh.x,mesh.y,mesh.z] = meshgrid(0: (1/numel(J_set)): J_set(length(J_set)), 0: (1/numel(D_set)): D_set(length(D_set)), 0: (1/numel(m_set)): m_set(length(m_set)));
Plot = slice(Valid_sol_counter_JDm, 1:(length(J_set)),1:(length(D_set)), 1:(length(m_set)));
colormap(mymap)
hold on
%plot3(critical_D(:), critical_J(:), critical_m(:), 'k.');
xlabel('J / U');
ylabel('\Delta / U');
zlabel('m');
set(Plot, 'EdgeColor','none', 'FaceColor','interp');
alpha(.1);
XT = get(gca, 'XTick');
set(gca, 'XTick', XT, 'XTicklabel', XT / numel(J_set) * J_set(length(J_set)));
YT = get(gca, 'YTick');
set(gca, 'YTick', YT, 'YTicklabel', YT / numel(D_set) * D_set(length(D_set)));
ZT = get(gca, 'ZTick');
set(gca, 'ZTick', ZT, 'ZTicklabel', ZT / numel(m_set) * m_set(length(m_set)));
saveas(gcf, fullfile(fname, 'JDm.eps'), 'epsc'); 
saveas(gcf, 'JDm.pdf'); 

