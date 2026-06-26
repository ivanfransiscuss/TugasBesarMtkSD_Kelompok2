% Menghapus hasil output, variabel, dan figur sebelumnya
clc;
clear;
close all;

%% ==========
%% Bagian 3
%% ==========

%% Gambar MRI menjadi matriks A
A = double(imread('MRI_MTK_SD.png'));
% Mengubah gambar ke grayscale
if ndims(A) == 3
   A = rgb2gray(uint8(A));
   A = double(A);
end
[m,n] = size(A);
fprintf('Ukuran matriks A = %d x %d\n',m,n);
%% SVD
[U,S,V] = svd(A);
%% Nilai-nilai k
k_values = [5 10 20 50 100];
%% Tabel hasil
fprintf('\n');
fprintf('===============================================================\n');
fprintf(' k\tCR\t\tRMSE\t\tPSNR(dB)\n');
fprintf('===============================================================\n');
for idx = 1:length(k_values)
   k = k_values(idx);
   %% Rekonstruksi rank-k
   Uk = U(:,1:k);
   Sk = S(1:k,1:k);
   Vk = V(:,1:k);
   Ak = Uk * Sk * Vk';
   %% Compression Ratio
   CR = (256^2) / (k*(256+256+1));
   %% RMSE
   rmse = sqrt(mean((A(:)-Ak(:)).^2));
   %% PSNR
   psnr_value = 20*log10(255/rmse);
   %% Mean, Variance, Covariance untuk SSIM
   mux = mean(A(:));
   muy = mean(Ak(:));
   sigx2 = var(A(:));
   sigy2 = var(Ak(:));
   C = cov(A(:),Ak(:));
   sigxy = C(1,2);
   %% Menampilkan hasil utama
   fprintf('%d\t%.4f\t\t%.4f\t\t%.4f\n', ...
           k,CR,rmse,psnr_value);
   %% Menampilkan parameter SSIM
   fprintf('\n===== Parameter SSIM untuk k = %d =====\n',k);
   fprintf('Mean Original (mux)       = %.6f\n',mux);
   fprintf('Mean Rekonstruksi (muy)   = %.6f\n',muy);
   fprintf('Variance Original         = %.6f\n',sigx2);
   fprintf('Variance Rekonstruksi     = %.6f\n',sigy2);
   fprintf('Covariance                = %.6f\n',sigxy);
   fprintf('========================================\n\n');
   %% Menampilkan citra hasil rekonstruksi
   figure;
   imshow(uint8(Ak));
   title(['Rekonstruksi SVD, k = ',num2str(k)]);
end
%% Menampilkan citra asli
figure;
imshow(uint8(A));
title('Citra MRI Asli');
