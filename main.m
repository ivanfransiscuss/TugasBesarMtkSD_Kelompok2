% Menghapus hasil output, variabel, dan figur sebelumnya
clc;
clear;
close all;

%% ==========
%% BAGIAN 1
%% ==========

% load mri
img = imread('MRI MTK SD.png');

% jika RGB ubah menjadi grayscale
if size(img,3) == 3
    img = rgb2gray(img);
end

% konversi ke double
A = double(img);

fprintf('Ukuran Matriks A = %d x %d\n', size(A,1), size(A,2));

% verifikasi Rank
r = rank(A);

fprintf('Rank(A) = %d\n', r);

% ambil submatriks 4x4
A4 = A(100:103,100:103);

disp('Submatriks 4x4 :');
disp(A4);

% determinan 
detA4 = det(A4);

fprintf('Determinan A4 = %.0f\n', detA4);

%% ==========
%% BAGIAN 2
%% ==========

% load mri
img = imread('MRI MTK SD.png');

if size(img,3) == 3
    img = rgb2gray(img);
end

A = double(img);

% ambil submatriks 3x3
B = A(95:97,96:98);

disp('Submatriks B = ');
disp(B);

% determinan Manual
detB = ...
    B(1,1)*(B(2,2)*B(3,3)-B(2,3)*B(3,2)) ...
    - B(1,2)*(B(2,1)*B(3,3)-B(2,3)*B(3,1)) ...
    + B(1,3)*(B(2,1)*B(3,2)-B(2,2)*B(3,1));

fprintf('Determinan B = %.0f\n', detB);

% matriks Kofaktor

C11 = det([69 71;73 76]);
C12 = -det([56 71;59 76]);
C13 = det([56 69;59 73]);

C21 = -det([61 70;73 76]);
C22 = det([52 70;59 76]);
C23 = -det([52 61;59 73]);

C31 = det([61 70;69 71]);
C32 = -det([52 70;56 71]);
C33 = det([52 61;56 69]);

C = [C11 C12 C13;
    C21 C22 C23;
    C31 C32 C33];

disp('Matriks Kofaktor = ');
disp(C);

% matriks Adjoin
adjB = C';

disp('Matriks Adjoin = ');
disp(adjB);

% invers Metode Adjoin

Binv = (1/detB)*adjB;

disp('Invers B = ');
disp(Binv);

% verifikasi

Itest = B*Binv;

disp('B * B^-1 = ');
disp(Itest);

% error numerik

err = norm(Itest-eye(3),'fro');

fprintf('Error Frobenius = %.4e\n',err);

% bandingkan dengan inv()

Binv_matlab = inv(B);

disp('Invers MATLAB = ');
disp(Binv_matlab);

selisih = norm(Binv-Binv_matlab,'fro');

fprintf('Selisih terhadap inv() = %.4e\n',selisih);

%% ==========
%% Bagian 3
%% ==========

%% Gambar MRI menjadi matriks A
A = double(imread('MRI MTK SD.png'));

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

   %% Konstanta stabilitas SSIM (Standar pengolahan citra)
   L = 255;
   C1 = (0.01 * L)^2;
   C2 = (0.03 * L)^2;

   %% Rumus SSIM 
   ssim_value = ((2 * mux * muy + C1) * (2 * sigxy + C2)) / ((mux^2 + muy^2 + C1) * (sigx2 + sigy2 + C2));

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
   fprintf('Nilai Akhir SSIM          = %.6f\n',ssim_value);
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

%% ==========
%% Bagian 4
%% =========

%% Membaca citra MRI
A = imread('MRI MTK SD.png');
%% Cek lagi jika punya 3 channel (RGB) ubah ke grayscale
if ndims(A) == 3
    A = rgb2gray(A);
end
A = double(A);

%% Singular Value Decomposition
[U,S,V] = svd(A);

%% Aproksimasi rank-10
k = 10;
A10 = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';
k = 10;
error_numerik = norm(A-A10,'fro');
sigma = diag(S);
error_teori = sqrt(sum(sigma(k+1:end).^2));
fprintf("Error teori   = %.6f\n",error_teori);
fprintf("Error numerik = %.6f\n",error_numerik); 


%% Error Frobenius hasil SVD
error_svd = norm(A - A10,'fro');
fprintf('Error Frobenius (SVD) = %.4f\n', error_svd);

%% Membuat matriks rank-10 pembanding
B = rand(size(A,1),k) * rand(k,size(A,2));
error_random = norm(A - B,'fro');
fprintf('Error Frobenius (Random Rank-10) = %.4f\n', error_random);

%% Kesimpulan otomatis
if error_svd < error_random
    fprintf('\nHasil SVD memberikan error yang lebih kecil.\n');
    fprintf('Teorema Eckart-Young terverifikasi.\n');
else
    fprintf('\nPerlu dilakukan pengujian ulang.\n');
end

%% Menampilkan gambar
figure;
subplot(1,2,1)
imshow(A, [])
title('MRI Asli');

subplot(1,2,2)
imshow(A10, [])
title('Rekonstruksi Rank-10');

%% ==========
%% BAGIAN 5
%% ==========

%% Membaca citra MRI
A = imread('MRI MTK SD.png');

% Jika RGB ubah ke grayscale
if size(A,3)==3
    A = rgb2gray(A);
end

A = double(A);

%% Membuat koordinat piksel
[x,y] = meshgrid(1:size(A,2),1:size(A,1));

%% Parameter noise
Amp = 20;
wx = 0.2;
wy = 0.2;

%% Model noise sinusoidal
f = Amp*sin(wx*x).*sin(wy*y);

%% Menambahkan noise ke MRI
A_noise = A + f;

%% Turunan parsial pertama terhadap x
dfdx = Amp*wx*cos(wx*x).*sin(wy*y);

%% Turunan kedua terhadap x
d2fdx2 = -Amp*(wx^2)*sin(wx*x).*sin(wy*y);

%% Turunan parsial pertama terhadap y
dfdy = Amp*wy*sin(wx*x).*cos(wy*y);

%% Magnitude gradien
grad = sqrt(dfdx.^2 + dfdy.^2);

%% Menampilkan MRI asli
figure;
imshow(uint8(A));
title('MRI Asli');

%% Menampilkan noise sinusoidal
figure;
imshow(f,[]);
title('Noise Sinusoidal');

%% Menampilkan MRI + Noise
figure;
imshow(A_noise,[]);
title('MRI + Noise');

%% Menampilkan turunan parsial df/dx
figure;
imshow(dfdx,[]);
title('Turunan Parsial df/dx');

%% Menampilkan turunan kedua d2f/dx2
figure;
imshow(d2fdx2,[]);
title('Turunan Kedua d2f/dx2');

%% Menampilkan magnitude gradien
figure;
imshow(grad,[]);
title('Magnitude Gradien');

%% ==========
%% Bagian 6
%% ==========

%% Membaca citra MRI
A = imread('MRI MTK SD.png');

if size(A,3)==3
    A = rgb2gray(A);
end

A = double(A);

%% Menentukan ROI
ROI = A(80:180,80:180);

%% Menampilkan ROI
figure;
imshow(uint8(ROI));
title('Region of Interest (ROI)');

%% Membuat koordinat
[x,y] = meshgrid(1:size(ROI,2),1:size(ROI,1));

%% Parameter noise
Amp = 20;
wx = 0.2;
wy = 0.2;

%% Fungsi noise
f = Amp*sin(wx*x).*sin(wy*y);

%% Visualisasi fungsi noise
figure;
imshow(f,[]);
title('Fungsi Noise f(x,y)');

%% METODE RIEMANN

IR = sum(f(:));

%% METODE TRAPEZOID

IT = trapz(trapz(f));

%% METODE SIMPSON

temp = zeros(size(f,1),1);

for i = 1:size(f,1)

    baris = f(i,:);

    n = length(baris);

    if mod(n-1,2)~=0
        baris = baris(1:end-1);
        n = n-1;
    end

    h = 1;

    hasil = baris(1) + baris(end);

    for j = 2:n-1
        if mod(j,2)==0
            hasil = hasil + 4*baris(j);
        else
            hasil = hasil + 2*baris(j);
        end
    end

    temp(i) = h/3 * hasil;

end

n = length(temp);

if mod(n-1,2)~=0
    temp = temp(1:end-1);
    n = n-1;
end

hasil = temp(1) + temp(end);

for i = 2:n-1
    if mod(i,2)==0
        hasil = hasil + 4*temp(i);
    else
        hasil = hasil + 2*temp(i);
    end
end

IS = (1/3)*hasil;

%% HASIL

fprintf('\n');
fprintf('=============================\n');
fprintf('HASIL INTEGRAL NUMERIK\n');
fprintf('=============================\n');
fprintf('Riemann    = %.4f\n',IR);
fprintf('Trapezoid  = %.4f\n',IT);
fprintf('Simpson    = %.4f\n',IS);
fprintf('=============================\n');

%% Analisis Konvergensi

selisih_RT = abs(IR-IT);
selisih_TS = abs(IT-IS);

fprintf('\n');
fprintf('Selisih Riemann-Trapezoid = %.4f\n',selisih_RT);
fprintf('Selisih Trapezoid-Simpson = %.4f\n',selisih_TS);

%% Grafik Perbandingan

nilai = [IR IT IS];

figure;
bar(nilai);

set(gca,'XTickLabel',{'Riemann','Trapezoid','Simpson'});

ylabel('Nilai Integral');

title('Perbandingan Metode Integral Numerik');

grid on;

%% ==========
%% Bagian 7
%% ==========

%% Membaca citra MRI
A = imread('MRI MTK SD.png');
if ndims(A) == 3
   A = rgb2gray(A);
end
A = double(A);

%% SVD
[U,S,V] = svd(A);

%% Norma Frobenius citra asli
normA = norm(A,'fro');

%% Rank maksimum matriks
r = rank(A);

%% Menghitung error relatif
errorRelatif = zeros(r,1);
for k = 1:r
   Ak = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';
   errorRelatif(k) = norm(A-Ak,'fro') / normA;
end

%% Menentukan nilai k
k10 = find(errorRelatif < 0.10,1,'first');
k5  = find(errorRelatif < 0.05,1,'first');
k1  = find(errorRelatif < 0.01,1,'first');

%% Menampilkan hasil
fprintf('=====================================\n');
fprintf('Hasil Analisis Konvergensi\n');
fprintf('=====================================\n\n');
fprintf('Error relatif < 10%% dicapai pada k = %d\n',k10);
fprintf('Error relatif < 5%%  dicapai pada k = %d\n',k5);
fprintf('Error relatif < 1%%  dicapai pada k = %d\n',k1);

%% Grafik Konvergensi
figure;
plot(1:r,errorRelatif,'LineWidth',2);
hold on;
yline(0.10,'--');
yline(0.05,'--');
yline(0.01,'--');
plot(k10,errorRelatif(k10),'o');
plot(k5,errorRelatif(k5),'o');
plot(k1,errorRelatif(k1),'o');
grid on;
xlabel('Nilai k');
ylabel('Error Relatif');
title('Konvergensi Error Relatif Aproksimasi SVD');
legend('Error Relatif','10%','5%','1%');
