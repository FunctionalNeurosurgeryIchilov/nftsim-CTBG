clear all; close all; 
configuration = 'abeysuria2015'; % e-erps-all-nodes  eirs-corticothalamic ei-cortical abeysuria2015 fit-braintrak-reproduce_00

nf_struct = run_nftsim_pwd(configuration);
nf.report(nf_struct);
trace = nf.extract(nf_struct, {'Propagator.1.phi'});
[data, longside_nodes, shortside_nodes] = nf.grid(nf_struct, 'Propagator.1.phi');

% nf.get_frequencies
% nf.rfft

nf.movie(nf_struct, 'Propagator.1.phi', 0);
nf.plot_timeseries(nf_struct,  {'Propagator.1.phi'}, {[1:256]}, false, false);
[f_spatial, P_spatial, Kx, Ky, Pkf, x, y, Prf] = nf.spatial_spectrum(nf_struct,'Propagator.1.phi');
figure;loglog(f_spatial,P_spatial);title('power spectrum similar to analytic');
[f, P] = nf.spectrum(nf_struct,'Propagator.1.phi');
figure;loglog(f,P);

fs = 2*max(f);
t = 0:(1/fs):2;
x = chirp(t,3,1,8,'quadratic');
cfreqs = linspace(1, 10, 100);
nf.wavelet_spectrogram(nf_struct,'Propagator.1.phi', x, fs, cfreqs, 3, 'plot');
