% JkD
fname = 'C:\von_Server\ETH\BSc Physics\7\Bachelorarbeit\plots';
U = 1;
k_set = linspace(0,10,10);
D_set = linspace(0,10,numel(k_set));
J_set = linspace(0,10,numel(k_set));
cyan = [0 1 1];
yellow = [1 1 0];


figure
syms k J D
fnc = sqrt(3)/2 *k - D - J;
fs = fimplicit3(fnc, [0 10 0 10 0 10]);
fs.EdgeColor = 'none';
fs.FaceColor = 'cyan';
fs.FaceAlpha = 0.5;
xlabel('\kappa / U');
ylabel('J / U');
zlabel('\Delta / U');
hold on
h1 = fill3([8.66, 0, 8.66], [8.66, 8.66, 0], [0, 0, 0], 'c');
h1.EdgeColor = 'none';
f1.FaceAlpha = 0.5;
hold on
h2 = fill3([0, 8.66, 0], [0, 0, 8.66], [0, 0, 0], 'c');
h2.EdgeColor = 'none';
h2.FaceAlpha = 0.5;
h3 = fill3([0 8.66 8.66], [8.66 8.66 8.66], [0 0 10], 'c');
h3.EdgeColor = 'none';
h3.FaceAlpha = 0.5;
h4 = fill3([8.66 0 0], [8.66 8.66 8.66], [10 10 0], 'c');
h4.EdgeColor = 'none';
h4.FaceAlpha = 0.5;
h5 = fill3([8.66 8.66 8.66], [0 8.66 8.66], [0 0 10], 'c');
h5.EdgeColor = 'none';
h5.FaceAlpha = 0.5;
h6 = fill3([8.66 8.66 8.66], [8.66 0 0], [10 0 10], 'c');
h6.EdgeColor = 'none';
h6.FaceAlpha = 0.5;


h7 = fill3([0 8.66 0], [0 0 0], [0 10 10], 'y');
h7.EdgeColor = 'none';
h7.FaceAlpha = 0.5;
h8 = fill3([0 0 0], [0 8.66 0], [0 10 10], 'y');
h8.EdgeColor = 'none';
h8.FaceAlpha = 0.5;
h9 = fill3([0 0 8.66], [0 8.66 0], [10 10 10], 'y');
h9.EdgeColor = 'none';
h9.FaceAlpha = 0.5;

h10 = fill3([0 8.66 8.66], [0 0 0], [0 0 10], 'c');
h10.EdgeColor = 'none';
h10.FaceAlpha = 0.5;
h11 = fill3([0 0 0], [0 8.66 8.66], [0 0 10], 'c');
h11.EdgeColor = 'none';
h11.FaceAlpha = 0.5;
%[mesh_k,mesh_J,mesh_D] = meshgrid(0: (1/numel(k_set)): k_set(length(k_set)), 0: (1/numel(J_set)): J_set(length(J_set)), 0: (1/numel(D_set)): D_set(length(D_set)));
saveas(gcf, fullfile(fname, 'kJD_condition.eps'), 'epsc'); 
saveas(gcf, 'kJD_condition.pdf'); 


%fill3([0 8 10], [10, 10, 0], [8, 0, 10], 'c');
    
%kD_vec = [20/sqrt(3), 0, 10];
%JD_vec = [0, ]

%{
    Matrix = zeros(numel(k_set), numel(k_set), numel(k_set));
    f = @(k, J, D) (sqrt(3)/2 *k - D - J);
    for Ji = 1: numel(J_set)
        J = J_set(Ji);
        for ki = 1:numel(k_set)
            k = k_set(ki);
            for Di = 1:numel(D_set)
                D = D_set(Di);
                if (double(subs(fnc)) + J < 0) && (double(subs(fnc)) + D < 0) && (double(subs(fnc)) - sqrt(3)/2*k > 0)
                    Matrix(k+1, J+1, D+1) = 1;
                else
                    Matrix(k+1, J+1, D+1) = 3;
                end
            end
        end
    end
    
    slice(Matrix1:(length(k_set)),1:(length(J_set)), 1:(length(D_set)));
    
  %}      



