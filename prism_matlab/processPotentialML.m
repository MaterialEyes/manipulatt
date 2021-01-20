 function [sTraining] = processPotentialML(emdSTEM)

% Colin Ophus - 2020 Dec
% Generate potential images for tilted samples for ML AI training.

% Inputs:
% emdSTEM - structure with simulation data
% zInds - indices for z direction integration of the potentials
qzSigma = 1/10;  % Band limiting of the potential array  (1/Angstroms)

% Plotting Inputs: (none of these variables affect the output)
flagPlotScatter = false;   % plot the potentials / qz as scatter plots
flagPlotImage = false;     % plot the potentials / qz as an image

% Outputs
% sTraining    - struct containing all training data:
%    dataPot   - [Nx Ny Nt] array, potential Fourier coefficients (Volts)
%    dataQz    - [Nx Ny Nt] array, out-of-plane displacements  (1/Ang)
%    dataMeas  - [Nx Ny Nt] array, probe images measured after scattering 
%                           through sample, for each of Nt total thicknesses.
%    dataProbe - [Nx Ny] Initial probe image (corresponding to vacuum image
%                        in experiments).  Both of these images normalized
%                        to an initial maximum value of the probe = 1.  
%                        (we should consider other normalizations too)
%    qx        - qx values (1/Ang)
%    qy        - qy values (1/Ang)
%    thickness - output thicknesses (Ang)

% input variables
numOutputPlanes = length(emdSTEM.output4Dthicknesses);
sliceThickness = emdSTEM.sliceThickness;
pixelSizeAA = 2 * emdSTEM.pixelSize;  % scaling *2 due to anti-aliasing

% init coordinates
N = [size(emdSTEM.pot,1)/2 size(emdSTEM.pot,2)/2 numOutputPlanes];
qx = reshape(fftshift(makeFourierCoords(N(1),pixelSizeAA(1))),[N(1) 1 1]);
qy = reshape(fftshift(makeFourierCoords(N(2),pixelSizeAA(2))),[1 N(2) 1]);

% initialize output struct
sTraining.dataPot  = zeros(N,'single');
sTraining.dataQz   = zeros(N,'single');
sTraining.dataMeas = zeros(N,'single');
sTraining.dataProbe = zeros(N(1:2),'single');
sTraining.qx = qx;
sTraining.qy = qy;
sTraining.thickness = emdSTEM.output4Dthicknesses;

% initial probe intensity / vacuum probe measurement
intProbe = fftshift(abs(emdSTEM.PsiInit(emdSTEM.xAA,emdSTEM.yAA)).^2);
intScale = 1 / max(intProbe(:));
sTraining.dataProbe(:) = intProbe * intScale;

% Output measurements - scale intensities to match probe
sTraining.dataMeas(:) = ...
    intScale * fftshift(fftshift( ...
    emdSTEM.output4D(:,:,1,1,:), 1), 2);


% initial windowing functions in x and y
wx = reshape(tukeywinMake(2*N(1)),[2*N(1) 1 1]);
wy = reshape(tukeywinMake(2*N(2)),[1 2*N(2) 1]);

% scaling parameter for outputs
scale = (1/4) / prod(N(1:2));


% main loop over all thicknesses
for a0 = 1:numOutputPlanes
    % subset of the potentials corresponding to this thickness
    indsRange = 1:emdSTEM.output4Dinds(a0);
    numPlanes = length(indsRange);
    
    
    % windowing function and filter over z dimension
    wz = reshape(tukeywinMake(numPlanes),[1 1 numPlanes]);
    qz = reshape(makeFourierCoords(numPlanes,sliceThickness),[1 1 numPlanes]);
    qzFilter = exp(-qz.^2/(2*qzSigma^2));

    % windowed FFT
    potFFT = abs( ...
        fftn(emdSTEM.pot(:,:,indsRange) ...
        .* (scale*wx.*wy.*wz))) .* qzFilter;
    
    % sum of peaks, qz estimate for all pixels
    potSum = fftshift(sum(potFFT(emdSTEM.xAA,emdSTEM.yAA,:),3));
    qzPot = fftshift(sum(potFFT(emdSTEM.xAA,emdSTEM.yAA,:) .* qz,3)) ./ potSum;
    
    % Output values
    sTraining.dataPot(:,:,a0) = potSum;
    sTraining.dataQz(:,:,a0) = qzPot;
    
    %     % testing for tungsten cell - mean inner potential ~37 Volt
    %     vx = (N(1)/2+1) + (-3:3);
    %     vy = (N(2)/2+1) + (-3:3);
    %     potTotal = sum(sum(potSum(vx,vy)));
    %     estimated_thickness = potTotal / 37;
    %     [sTraining.thickness(a0) estimated_thickness]
end


 if  (flagPlotImage == true) || (flagPlotScatter == true)
     % plotting colormaps
     % colormap - blue-white-red
     c = linspace(0,1,64)';
     c0 = zeros(64,1);
     c1 = ones(64,1);
     cmap = [ ...
         [c0 0.2+0.4*c 1-0.4*flipud(c)];
         [c 0.6+0.4*c c1];
         [c1 1-c 1-c];
         [1-0.6*c c0 c0];
         ];
     x = linspace(0,1,128)';
     y = 5*x.^4 - 4*x.^5;
     dc = [y; flipud(y)];
     w = 0.8;
     cmap(:) =  cmap.*(1 - dc*w) + dc*w;
 end
 
 indPlot = 21

 
 if flagPlotImage == true
     % images
     
     %     indPlot = numOutputPlanes;
    % indPlot = 9;
     
     powerScale = 0.5;
     intRangeProbes = [0 0.1];
     intRangePot = [0.2 2];
     intRange_qz = [-1 1]*0.2;  % mrads
     
     % Generate colour image to show pot values and qz
     Iv = (sTraining.dataPot(:,:,indPlot) ...
         / sTraining.thickness(indPlot)) .^ powerScale;
     Iv(:) = (Iv - intRangePot(1)) / (intRangePot(2) - intRangePot(1));
     Iv(:) = min(max(Iv,0),1);
     Iqz = (sTraining.dataQz(:,:,a0) - intRange_qz(1)) ...
         / (intRange_qz(2) - intRange_qz(1));
     Iqz(:) = min(max(Iqz,0),1);
     Irgb = ind2rgb(round(255*Iqz)+1,cmap);
     Irgb = Irgb .* Iv;
     
     % STEM probes (initial and measured)
     fig1 = figure(101)
     clf
     imagesc([sTraining.dataProbe ...
         sTraining.dataMeas(:,:,indPlot)] .^ powerScale)
     axis equal off
     colormap(gray(256))
     colorbar
     caxis(intRangeProbes)
     
     saveas(fig1,'trainset_output1','png') 
    
     % potential and qz
     fig2 = figure(102)
     clf
     imagesc(Irgb)
     axis equal off
     colormap(gray(256))
     colorbar
     saveas(fig2,'trainset_output2','png')
     %     caxis(intRangePot)
     
 end
 
 
 if flagPlotScatter == true
     % overlaid scatter plot
    % indPlot = 10;
     powerScale = 0.5;
     intRangeProbes = [0 0.1];
     potMinPlot = 0.25;
     intRange_qz = [-1 1]*0.2;  % mrads
     
     % find peaks to plot
     im = sTraining.dataPot(:,:,indPlot) / sTraining.thickness(indPlot);
     p = im > circshift(im,[-1 -1]) & ...
         im > circshift(im,[ 0 -1]) & ...
         im > circshift(im,[ 1 -1]) & ...
         im > circshift(im,[-1  0]) & ...
         im > circshift(im,[ 1  0]) & ...
         im > circshift(im,[-1  1]) & ...
         im > circshift(im,[ 0  1]) & ...
         im > circshift(im,[ 1  1]) & ...
         im > potMinPlot;
     [xp,yp] = find(p);
     qxp = qx(xp);
     qyp = qy(yp);
     Iqz = sTraining.dataQz(:,:,indPlot);
     qzp = Iqz(p);
     potSig = im(p);
         
     Iv = sTraining.dataMeas(:,:,indPlot) .^ powerScale;
     Iv(:) = (Iv - intRangeProbes(1)) / (intRangeProbes(2) - intRangeProbes(1));
     Iv(:) = min(max(Iv,0),1);
     
     mSizeScale = 60;
     
     fig3 = figure(103)
     clf
     imagesc( ...
         repmat(Iv,[1 1 3]),...
         'xdata',qy(:),...
         'ydata',qx(:))
     hold on
     scatter( ...
         qyp,...
         qxp,...
         potSig * mSizeScale,...
         qzp,...
         'filled',...
         'marker','o',...
         'linewidth',1,...
         'markeredgecolor',[0 0 0])
     
     hold off
     axis equal off
     colormap(cmap)
     %     colorbar
     caxis(intRange_qz)
     set(gca,'position',[0 0 1 1])
     saveas(fig3,'trainset_output3','png')
   
end



end
