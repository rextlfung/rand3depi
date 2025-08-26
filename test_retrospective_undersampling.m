%% Short script for testing retrospective undersampling
load("/mnt/storage/rexfung/20250819retroball/gt_estimate.mat")
load("/mnt/storage/rexfung/20250819retroball/smaps.mat")

%% Crop FOV in z (to be able to run recon on this machine)
gt_estimate = gt_estimate(:,:,25:(108-24));
[Nx, Ny, Nz] = size(gt_estimate);
smaps_raw = smaps_raw(:,:,25:(108-24),:);

%% Create artifial time series
Nt = 50;
img = repmat(gt_estimate, [1, 1, 1, Nt]);

%% Generate random sampling mask
Ry = sqrt(3); Rz = sqrt(3); R = [Ry Rz];
caipi_z = 2;
acs_y = 0; acs_z = 0; acs = [acs_y acs_z];
weights_y = normpdf(1:Ny, mean(1:Ny), Ny/6);
weights_z = normpdf(1:Nz, mean(1:Nz), Nz/6);
max_ky_step = round(Ny/16);

omegas = zeros(Ny, Nz, Nt);
for t = 1:Nt
    rng(t);
    [omegas(:,:,t), ~, ~]= randsamp2dcaipi([Ny, Nz], R, acs, weights_y, weights_z, max_ky_step, caipi_z);
end

%% Generate multicoil images
img_mc = smaps_raw .* reshape(img,Nx,Ny,Nz,1,Nt);

%% 3D FFT to get k-space
ksp_mc = zeros(size(img_mc));
for t = 1:Nt
    ksp_mc(:,:,:,:,t) = toppe.utils.ft3(img_mc(:,:,:,:,t));
end
clear img_mc;

%% Undersample
ksp_us_mc = ksp_mc .* reshape(omegas,1,Ny,Nz,1,Nt);
clear ksp_mc;
save("/mnt/storage/rexfung/20250819retroball/ksp_us_6x.mat", "ksp_us_mc", "-v7.3");

%% 3D IFFT to get multicoil aliased images
img_mc = zeros(size(ksp_us_mc));
for t = 1:Nt
    img_mc(:,:,:,:,t) = toppe.utils.ift3(ksp_us_mc(:,:,:,:,t));
end
clear ksp_us_mc;

%% Sensitivity map weighted coil combination
img_final = squeeze(sum(conj(smaps_raw) .* img_mc, 4) ./ (sum(abs2(smaps_raw), 4) + eps));
clear img_mc;

%% Scrollable plot
interactive4D(abs(permute(img_final, [2 1 3 4])));

%% tSNR maps
tSNR_map = mean(abs(img_final), 4) ./ (std(abs(img_final), 0, 4) + eps);
figure; im(tSNR_map); colormap hot; colorbar;