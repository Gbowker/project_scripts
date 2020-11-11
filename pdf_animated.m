% define alpha...
csHex = crystalSymmetry('6/mmm', [3 3 4.7], 'X||b', 'Y||a*', 'Z||c*', 'mineral', 'Ti-Hex', 'color', [0.53 0.81 0.98]);
       %crystalSymmetry('6/mmm', [3 3 4.7], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Ti-Hex', 'color', [0.53 0.81 0.98]);
% define beta...
csCubic = crystalSymmetry('m-3m', [3.2 3.2 3.2], 'mineral', 'Titanium cubic', 'color', [0.56 0.74 0.56]);


% Need to change cs if alpha or beta:
cs = csHex;

% Define path to ori_ files:
path_to_incs = "../DAMASK_results/mixed_1024grains/alpha/";
mkdir ../DAMASK_results/mixed_1024grains/alpha/all_10-10_poleplots/
mkdir ../DAMASK_results/mixed_1024grains/alpha/all_odfs/

%% loop over increments...

for inc = 0.0:200.0
    
    increment = string(inc);
    disp(inc)
    strain = string(round(inc/2.0, 2));
    
    % startupDefine path to quaternion.txt files...
    path = strcat(path_to_incs, 'ori_inc', increment, '.txt');
    
    % Read the quaternions from the data file
    fid = fopen(path);
    data = textscan(fid, '%f%f%f%f', 'HeaderLines', 1, 'CollectOutput', 1);
    data = data{:};
    fid = fclose(fid);
    q = quaternion(transpose(data));
    
    % Estiamte an ODF from the orientations
    ori = orientation(q, cs);
    % calcdensity method
    odf = calcDensity(ori,'kernel',deLaValleePoussinKernel,'halfwidth',5*degree);
    % elliot et als help
%     psi = calcKernel(ori,'Method','KLCV')
%     odf = calcKernelODF(ori,'kernel',psi);
    
    % plot pole figures
    f1 = figure('visible','off');
    plotPDF(odf, Miller(1, 0, -1, 0, cs),'antipodal', 'contourf', 0:0.1:10, 'minmax', 'colorRange',[0,2.0]); % nicks consistent colobar!
    annotate(xvector,'label',{'RD'},'BackgroundColor','w','MarkerSize',17,'FontSize', 26)
    annotate(yvector,'label',{'TD'},'BackgroundColor','w','MarkerSize',17,'FontSize', 26)
    annotate(zvector,'label',{'ND'},'BackgroundColor','w','MarkerSize',17,'FontSize', 26)
    mtexColorbar('location','southoutside','FontSize', 26)
    title(strcat('(10-10) % reduction = ', strain, '%'))
    saveas(gcf, strcat(path_to_incs,'all_10-10_poleplots/', increment,'_inc', '_poleplot.png'))
    close(f1)
    
    % plot odf phi2 slices
%     odf.SS=specimenSymmetry('orthorhombic');
%     f2 = figure('visible','off');
%     plotSection(odf,'contourf','phi2',45*degree);
%     mtexColorbar('location','southoutside','FontSize', 26)
%     saveas(gcf, strcat(path_to_incs,'all_odfs/inc_', increment, '_odfplot.png'))
%     close(f2)
end
