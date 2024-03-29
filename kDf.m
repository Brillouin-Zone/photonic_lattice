%%% spinodals in a  3D (k, D, f) plot
% MFT: J = 1; U = const

fname = 'C:\von_Server\ETH\BSc Physics\7\Bachelorarbeit\plots';
U = 1;
k_set = linspace(0,2,100);%[0.1 : (2/100): 2];
f_set = linspace(0,10,numel(k_set));%[0.1 : (10 / numel(k_set)) : 10];
D_set = linspace(0,10,numel(k_set));%[0.1 : (10 / numel(k_set)) : 10];
Valid_sol_counter = zeros(numel(k_set), numel(k_set), numel(k_set)); 
critical_k = zeros(1, numel(k_set));
critical_f = zeros(1, numel(k_set));
critical_D = zeros(1, numel(k_set));

for fi = 1: numel(f_set)
    f = f_set(fi);
    for Di = 1:numel(D_set)
        D = D_set(Di);
        for ki = 1:numel(k_set)
            k = k_set(ki);
            P = [1, -(2.*D./U), ((D.^2 + k.^2./4)/U.^2), -(f./U)^2];
            N0 = roots(P);
            N_set1(ki, Di, fi) = N0(1);
            N_set2(ki, Di, fi) = N0(2);
            N_set3(ki, Di, fi) = N0(3);
            
            if real(N_set1(ki, Di, fi)) > 0 && (imag(N_set1(ki, Di, fi)) == 0)
                Valid_sol_counter(ki, Di, fi) = Valid_sol_counter(ki, Di, fi) +1;
            end
            if real(N_set2(ki, Di, fi)) > 0 && (imag(N_set2(ki, Di, fi)) == 0)
                Valid_sol_counter(ki, Di, fi) = Valid_sol_counter(ki, Di, fi) +1;
            end
            if real(N_set3(ki, Di, fi)) > 0 && (imag(N_set3(ki, Di, fi)) == 0)
                Valid_sol_counter(ki, Di, fi) = Valid_sol_counter(ki, Di, fi) +1;
            end
        end
    end
end
%Valid_sol_counter = rot90(fliplr(Valid_sol_counter));

b = 0;
for ki = 1:numel(k_set) % kappa ebenen
    for Di = 1:numel(D_set) % 
        for fi = 1:numel(f_set)
            if Valid_sol_counter(ki, Di, fi) == 3
                critical_k(ki) = ki;
                critical_f(ki) = fi;
                critical_D(ki) = Di;
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
%Valid_sol_counter = rot90(fliplr(Valid_sol_counter));


figure
[mesh.x,mesh.y,mesh.z] = meshgrid(0: (1/numel(k_set)): k_set(length(k_set)), 0: (1/numel(D_set)): D_set(length(D_set)), 0: (1/numel(f_set)): f_set(length(f_set)));
Plot = slice(Valid_sol_counter, 1:(length(k_set)),1:(length(D_set)), 1:(length(f_set)));
hold on
plot3(critical_D(:), critical_k(:), critical_f(:), 'k.');
ylabel('\kappa / U');
xlabel('\Delta / U');
zlabel('f / U');
%title('kDf');
set(Plot, 'EdgeColor','none', 'FaceColor','interp');
alpha(.1);
XT = get(gca, 'XTick');
set(gca, 'XTick', XT, 'XTicklabel', XT / numel(D_set) * D_set(length(D_set)));
YT = get(gca, 'YTick');
set(gca, 'YTick', YT, 'YTicklabel', YT / numel(k_set) * k_set(length(k_set)));
ZT = get(gca, 'ZTick');
set(gca, 'ZTick', ZT, 'ZTicklabel', ZT / numel(f_set) * f_set(length(f_set)));
saveas(gcf, fullfile(fname, 'kDf.eps'), 'epsc'); 
saveas(gcf, 'kDf.pdf'); 



