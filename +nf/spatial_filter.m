%% Given a grid of voltages, applies a spatial filter that imitates volume conduction
%
% ARGUMENTS:
%        obj -- nf object | a 3D matrix of values
%        p -- trace to use | sampling rate
%        kmax -- is the number of k-values to use. This should match with whatever kmax is set to in the analytic summation
%
% OUTPUT:
%        data -- filtered grid of voltages
%        Kx -- vector with angular wavenumbers in [rad/m]
%        Ky -- vector with angular wavenumbers in [rad/m]
%        x -- array of spatial coordinates along x [m]
%        y -- array of spatial coordinates along y [m]
%
% REQUIRES:
%        nf.get_frequencies() 
%        nf.grid()  
%
% REFERENCES:
%
% AUTHOR:
%     Daniel Polyakov (2023-06-18).
%
% USAGE:
%{
    %
    data = nf.spatial_filter(nf, p, kmax)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data, Kx, Ky, x, y] = spatial_filter(obj, p, kmax)

    if nargin < 3 || isempty(kmax)
        kmax = []; % Use all k-values by default
    end

    % work out the sampling rate in pixels per metre
    if isstruct(obj)
        if nargin < 2 || isempty(p)
            p = 'propagator.1.phi'; % Try the phi propagator first
        end
        data = nf.grid(obj, p);
        fs = 1.0 / obj.deltat;
    else
        data = obj;
        fs = p;
    end

    if any(isnan(data(:)))
        data = data(:, :, 1:(end - 1));
        if any(isnan(data(:)))
            error('NaNs are present in the data set')
        end
    end

    if mod(size(data, 1), 2) || mod(size(data, 2), 2)
        error('In order to have a zero frequency component, you need an even number of grid edge nodes');
    end


    % Calculate the Fourier f and k values
    Lx = 0.5; % linear cortex dimension (m)
    Ly = 0.5;
    [~, Kx, Ky, x, y, ~, ~, ~] = nf.get_frequencies(data, fs, Lx, Ly);

    k2 = Kx.^2 + Ky.^2; % Matrix of k-squared values
    if isempty(kmax)
        k_mask = ones(size(k2));
    else
        [center_x, center_y] = find(k2 == 0); % Get the centre entry
        [a, b] = meshgrid(1:size(k2, 1), 1:size(k2, 2));
        k_mask = abs(a - center_x) <= kmax & abs(b - center_y) <= kmax;
    end
    
    % Calculate the value of k^2 at each grid node for spatial filtering
    k0 = 10; % spatial filter constant (m^-1)
    k_filter = exp(-k2 ./ k0^2);

%     data = bsxfun(@times, data, k_mask);
%     data = bsxfun(@times, data, sqrt(k_filter));

    % First, get the 3D FFT and apply volume conduction
    P = fftshift(fftn(data));
    P = bsxfun(@times, P, k_mask);
    P = bsxfun(@times, P, sqrt(k_filter));
    data = ifftn(ifftshift(P));

%     Pkf = Pkf + Pkf_new;
%     P = squeeze(sum(sum(Pkf, 1), 2)) * dk * dk;
