D = 3; N0 = 1000; k = 10; p = 0.01;
A = small_world_manhattan(N0, k, D, p); N = size(A, 1);
[ss, nn] = count_distances(A);

disp('Using model c(s) ∝ s^{D-1}:');
[DVec, sMaxVec, D2Mat, codeCell, descriptionCell, DMat, objFuncMat] = est_corr_dim_3(ss, nn);
DCell = arrayfun(@(x) x, DVec, 'Uni', false);
disp([DCell(1:3), descriptionCell(1:3)']);

disp('Using model C(s) ∝ s^D:');
[DVec, sMaxVec, D2Mat, codeCell, descriptionCell, DMat, objFuncMat] = est_corr_dim_4(ss, nn);
DCell = arrayfun(@(x) x, DVec, 'Uni', false);
disp([DCell, descriptionCell']);