% kJm mit m=1+(2D)/U
 % note: for single cavity => cavity array: use the substitution D => D + J
    % (x,y,z) = (J,k,m)
fname = 'C:\von_Server\ETH\BSc Physics\7\Bachelorarbeit\plots';
U = 1;
f = 1;
k_set = linspace(0,2,100);
D_set = linspace(0,10,numel(k_set));
J_set = linspace(0,10,numel(k_set));
DJ_set = (D_set(:) + J_set(:)).';
m_set = 1+ 2*DJ_set ./U;
    % D = [(m-1) * U / 2 ] - J 
            % or:  D + J = [(m-1) * U / 2 ]
Valid_sol_counter_kJm = zeros(numel(k_set), numel(k_set), numel(k_set)); 
%critical_k = zeros(1, numel(k_set));
%critical_J = zeros(1, numel(k_set));
%critical_m = zeros(1, numel(k_set));

for Ji = 1: numel(J_set)
    J = J_set(Ji);
    for ki = 1:numel(k_set)
        k = k_set(ki);
        for mi = 1:numel(m_set)
            %D + J = (U/2) * (m_set(mi)-1);
            m = m_set(mi);
            P_kJm = [1, -(2.* ((m-1)*U/2 +J) ./U), (( ((m-1)*U/2 +J).^2 + k.^2./4)/U.^2), -(f./U)^2]; % THIS ONE
            N0_kJm = roots(P_kJm);
            N_set1_kJm(Ji, ki, mi) = N0_kJm(1);
            N_set2_kJm(Ji, ki, mi) = N0_kJm(2);
            N_set3_kJm(Ji, ki, mi) = N0_kJm(3);
            
            if real(N_set1_kJm(Ji, ki, mi)) > 0 && (imag(N_set1_kJm(Ji, ki, mi)) == 0)
                Valid_sol_counter_kJm(Ji, ki, mi) = Valid_sol_counter_kJm(Ji, ki, mi) +1;
            end
            if real(N_set2_kJm(Ji, ki, mi)) > 0 && (imag(N_set2_kJm(Ji, ki, mi)) == 0)
                Valid_sol_counter_kJm(Ji, ki, mi) = Valid_sol_counter_kJm(Ji, ki, mi) +1;
            end
            if real(N_set3_kJm(Ji, ki, mi)) > 0 && (imag(N_set3_kJm(Ji, ki, mi)) == 0)
                Valid_sol_counter_kJm(Ji, ki, mi) = Valid_sol_counter_kJm(Ji, ki, mi) +1;
            end
        end
    end
end
%Valid_sol_counter = rot90(fliplr(Valid_sol_counter));

% surface Jm as boundary of bistable==smooth / stable==sharp region are the
% critical points
%{ 
b = 0;
for Ji = 1:numel(J_set) % kappa ebenen
    for ki = 1:numel(k_set) % 
        for mi = 1:numel(m_set)
            if Valid_sol_counter_kJm(Ji, ki, mi) == 3
                critical_J(Ji) = Ji;
                critical_k(Ji) = ki;
                critical_m(Ji) = mi;
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


mymap = [0 1 1
    1 1 0]; % cyan, yellow
figure
[mesh.x,mesh.y,mesh.z] = meshgrid(0: (1/numel(J_set)): J_set(length(J_set)), 0: (1/numel(k_set)): k_set(length(k_set)), 0: (1/numel(k_set)): m_set(length(m_set)));
Plot = slice(Valid_sol_counter_kJm, 1:(length(J_set)),1:(length(k_set)), 1:(length(m_set)));
colormap(mymap)
hold on
%plot3(critical_k(:), critical_J(:), critical_m(:), 'k.');
ylabel('\kappa / U');
xlabel('J / U');
zlabel('m = 1+ 2\Delta / U');
set(Plot, 'EdgeColor','none', 'FaceColor','interp');
alpha(.1);
XT = get(gca, 'XTick');
set(gca, 'XTick', XT, 'XTicklabel', XT / numel(J_set) * J_set(length(J_set)));
YT = get(gca, 'YTick');
set(gca, 'YTick', YT, 'YTicklabel', YT / numel(k_set) * k_set(length(k_set)));
ZT = get(gca, 'ZTick');
set(gca, 'ZTick', ZT, 'ZTicklabel', ZT / numel(m_set) * m_set(length(m_set)));
hax = gca;
hax.YTickLabel = flipud(hax.YTickLabel);
saveas(gcf, fullfile(fname, 'kJm.eps'), 'epsc'); 
saveas(gcf, 'kJm.pdf'); 
