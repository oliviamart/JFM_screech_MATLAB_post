%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get eigenvalues in region of the spectrum of interest
% and sort them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Olivia Martin
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nRnT = 70*70; % number of points Nr*Nt
num_eigs_tot = 200; % number of eigenvalues computed
x_val = 4;
freq_vec = [0.00:0.0025:0.35];
input_dir = '../../';
out_dir = '../';

U_inf = 0.007;
M = 1.78;
mu = 1;
h = 0.010; b = 0.040; 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP OVER FREQUENCIES AND LOAD IN DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eigvals = zeros(num_eigs_tot, length(freq_vec));
eigvect = cell(1, length(freq_vec));

for i = 1:length(freq_vec)

    V_vec = NaN(nRnT,num_eigs_tot);
    eigvals_vec = zeros(num_eigs_tot, 1);

    % LOAD DATA 
    filename = strcat(input_dir, 'results_x', num2str(x_val), '/Nr70_Nt70_f', ...
        num2str(freq_vec(i)), '_mu1_x', num2str(x_val), '.mat');
    if isfile(filename)
        load(filename);
        numE1 = size(D,1);
        D1 = diag(D);
        V1 = V(1:nRnT,:);
        eigvals_vec(1:numE1) = D1;
        V_vec(:,1:numE1) = V1;
    end

    eigvals(:,i) = eigvals_vec;
    eigvect{i} = V_vec;
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORT EIGENVALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eigvals_tracked] = sort_eigenvalues(eigvals, eigvect, freq_vec, U_inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FURTHER SORT EIGENVALUES BASED ON WHEN THEY APPEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eigvals_tracked2 = eigvals_tracked;
eigvals_sorted = NaN(size(eigvals_tracked)) + 1i*NaN;
count = 1;
for f = 1:size(eigvals_tracked2,2)
    [idx1] = find(~isnan(eigvals_tracked2(:,f)));
    [~,idx] = sort(-1*real(eigvals_tracked2(idx1,f)));
    for i = 1:length(idx)
        eigvals_sorted(count,:) = eigvals_tracked2(idx1(idx(i)),:);
        count = count + 1;
        eigvals_tracked2(idx1(idx(i)),:) = NaN;
    end
end

%% SORT EIGENFUNCTIONS
numF = size(eigvals_sorted,2);
numMode = size(eigvals_sorted,1);
eigvect_sorted = cell(1, numF); 

for j = 1:numF
    eigvect_sorted{j} = NaN(nRnT, numMode);
    for i = 1:numMode        
        eig = eigvals_sorted(i, j);
        if ~isnan(eig)  
            idx = find(abs(eigvals(:,j) - eig) == min(abs(eigvals(:,j) - eig)));
            eigvect_sorted{j}(:,i) = eigvect{j}(:,idx);
        end
    end
end

%% SAVE EIGENVALUES AND EIGENFUNCTIONS

filename = strcat(input_dir, 'results_x', num2str(x_val), '/Nr70_Nt70_f', ...
        num2str(freq_vec(1)), '_mu', num2str(mu), '_x', num2str(x_val), '.mat');
R = load(filename).R;
Th = load(filename).Th;
D0 = load(filename).D0;
T0 = load(filename).T0;

St = get_St(freq_vec, M, h, b);
filename_out = strcat(out_dir, 'eigs_sorted_Mj_', num2str(M), '_xVal_', num2str(x_val), ...
    '_mu', num2str(mu), '.mat');

save(filename_out, 'R', 'Th', 'D0', 'T0', 'eigvals_sorted', ...
    'eigvect_sorted', 'St', 'freq_vec');

%% PLOT

% figure
% hold on
% for i = 1:size(eigvals_sorted,1)
%     scatter(real(eigvals_sorted(i,:)), freq_vec)
% end
% xlim([-3 3])
% ylim([0 0.3])
% 
% figure; 
% p = D0 * T0 * eigvect_sorted{3}(:,1);
% p = reshape_eig(p, 70, 70);
% 
% pcolor(R.*cos(Th), R.*sin(Th), abs(p))
% shading interp
