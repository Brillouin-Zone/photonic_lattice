%DJf

    % (x,y,z) := (J, D, f)
fname = 'C:\von_Server\ETH\BSc Physics\7\Bachelorarbeit\plots';
J_set = linspace(0, 10, 100);
f_set = linspace(0,10,numel(J_set));
D_set = linspace(0,10,numel(J_set));
U = 1;
k = 1;
Valid_sol_counter_JDf = zeros(numel(J_set), numel(J_set), numel(J_set)); 
%critical_f = zeros(1, numel(J_set));
%critical_J = zeros(1, numel(J_set));
%critical_D = zeros(1, numel(J_set));

for Ji = 1:numel(J_set)
    J = J_set(Ji);
    for Di = 1:numel(J_set)
        D = D_set(Di);
        for fi = 1:numel(J_set)
            f = f_set(fi);
            P_JDf = [1, -(2.*(D+J)./U), (((D+J).^2 + k.^2./4)/U.^2), -(f./U)^2];
            N0_JDf = roots(P_JDf);
            N_set1_JDf(Ji, Di, fi) = N0_JDf(1);
            N_set2_JDf(Ji, Di, fi) = N0_JDf(2);
            N_set3_JDf(Ji, Di, fi) = N0_JDf(3);
            
            if real(N_set1_JDf(Ji, Di, fi)) > 0 && (imag(N_set1_JDf(Ji, Di, fi)) == 0)
                Valid_sol_counter_JDf(Ji, Di, fi) = Valid_sol_counter_JDf(Ji, Di, fi) +1;
            end
            if real(N_set2_JDf(Ji, Di, fi)) > 0 && (imag(N_set2_JDf(Ji, Di, fi)) == 0)
                Valid_sol_counter_JDf(Ji, Di, fi) = Valid_sol_counter_JDf(Ji, Di, fi) +1;
            end
            if real(N_set3_JDf(Ji, Di, fi)) > 0 && (imag(N_set3_JDf(Ji, Di, fi)) == 0)
                %Valid_sol_counter_JDf(Di, Ji, fi) = Valid_sol_counter_JDf(Di, Ji, fi) +1;
                Valid_sol_counter_JDf(Ji, Di, fi) = Valid_sol_counter_JDf(Ji, Di, fi) +1;
            end
        end
    end
end

% surface between bistable==smooth and stable==sharp region consists of the critical points
%{
b = 0;
for Di = 1:numel(J_set) % kappa ebenen
    for Ji = 1:numel(J_set) % 
        for fi = 1:numel(f_set)
            if Valid_sol_counter_JDf(Di, Ji, fi) == 3
                critical_J(Ji) = Ji;
                critical_f(Ji) = fi;
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
%Valid_sol_counter_JDf = rot90(fliplr(Valid_sol_counter_JDf));

figure
[mesh.x,mesh.y,mesh.z] = meshgrid(0: (1/numel(J_set)): J_set(length(J_set)), 0: (1/numel(J_set)): D_set(length(D_set)), 0: (1/numel(J_set)): f_set(length(f_set)));
Plot = slice(Valid_sol_counter_JDf, 1:(length(J_set)),1:(length(D_set)), 1:(length(f_set)));
hold on
%plot3(critical_D(:), critical_J(:), critical_f(:), 'k.');
xlabel('J / U');
ylabel('\Delta / U');
zlabel('f / U');
set(Plot, 'EdgeColor','none', 'FaceColor','interp');
alpha(.1);
XT = get(gca, 'XTick');
set(gca, 'XTick', XT, 'XTicklabel', XT / numel(J_set) * J_set(length(J_set)));
YT = get(gca, 'YTick');
set(gca, 'YTick', YT, 'YTicklabel', YT / numel(D_set) * D_set(length(D_set)));
ZT = get(gca, 'ZTick');
set(gca, 'ZTick', ZT, 'ZTicklabel', ZT / numel(f_set) * f_set(length(f_set)));
saveas(gcf, fullfile(fname, 'JDf.eps'), 'epsc'); 
saveas(gcf, 'JDf.pdf'); 















