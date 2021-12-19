% runs NFTsim from every directory
function nf_struct = run_nftsim_pwd(configuration)

nftsim_path = '..\nft\nftsim';
in_fn = ['.\configs\' configuration '.conf'];
out_fn = ['.\output\' configuration '.output'];
addpath(nftsim_path);

%%%%%%%%%%%%%%%%%%%%%%%%%

work_dir = pwd;
cd(nftsim_path);
disp(['NFTsim data generation: ' configuration]);
[status, cmdout] = system(['.\bin\nftsim.exe -i ' in_fn ' -o ' out_fn]);
if status ~= 0
    error(cmdout);
end
% if exist(out_fn,'file') ~=2
%     % nf_struct = nf.run(in_fn, false, [nftsim_path '\bin\nftsim.exe']);
%     [status, cmdout] = system(['.\bin\nftsim.exe -i ' in_fn ' -o ' out_fn]);
%     if status ~= 0
%         delete(out_fn);
%         error(cmdout);
%     end
% end
nf_struct = nf.read(out_fn);
nf_struct.conf_file = configuration;
cd(work_dir);
