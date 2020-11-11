%% Define MTEX preferences...
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outofPlane');
setMTEXpref('figSize','large');
setMTEXpref('FontSize', 20);

setMTEXpref('defaultColorMap', 'Viridis');

%% === define alpha and beta: === %

% define alpha...
csHex = crystalSymmetry('6/mmm', [3 3 4.7], 'X||b', 'Y||a*', 'Z||c*', 'mineral', 'Ti-Hex', 'color', [0.53 0.81 0.98]);
% define beta...
csCubic = crystalSymmetry('m-3m', [3.2 3.2 3.2], 'mineral', 'Titanium cubic', 'color', [0.56 0.74 0.56]);

%% Define path to quaternion.txt files...
alpha_single_oris = strcat('../DAMASK_results/alpha_1024grains/oris/ori_inc47.txt');
alpha_mixed_oris = strcat('../DAMASK_results/mixed_1024grains/alpha/ori_inc47.txt');

beta_single_oris = strcat('../DAMASK_results/beta_1024grains/oris/ori_inc_49.txt');
beta_mixed_oris = strcat('../DAMASK_results/mixed_1024grains/beta/ori_inc_49.txt');

%% Read alpha_singlephase data
alpha_file = fopen(alpha_single_oris);
alpha_single_data = textscan(alpha_file, '%f%f%f%f', 'HeaderLines', 1, 'CollectOutput', 1);
alpha_single_data = alpha_single_data{:};
alpha_file = fclose(alpha_file);
q_alpha_single = quaternion(transpose(alpha_single_data));
% Estimate alpha_single ODF
alpha_single_ori = orientation(q_alpha_single, csHex);
alpha_single_odf = calcDensity(alpha_single_ori,'kernel',deLaValleePoussinKernel,'halfwidth',5*degree);

%% Read alpha_mixedphase data
alpha_mixed_file = fopen(alpha_mixed_oris);
alpha_mixed_data = textscan(alpha_mixed_file, '%f%f%f%f', 'HeaderLines', 1, 'CollectOutput', 1);
alpha_mixed_data = alpha_mixed_data{:};
alpha_mixed_file = fclose(alpha_mixed_file);
q_alpha_mixed = quaternion(transpose(alpha_mixed_data));
% Estiamte alpha_single ODF
alpha_mixed_ori = orientation(q_alpha_mixed, csHex);
alpha_mixed_odf = calcDensity(alpha_mixed_ori,'kernel',deLaValleePoussinKernel,'halfwidth',5*degree);

%% Read beta_singlephase data
beta_file = fopen(beta_single_oris);
beta_single_data = textscan(beta_file, '%f%f%f%f', 'HeaderLines', 1, 'CollectOutput', 1);
beta_single_data = beta_single_data{:};
beta_file = fclose(beta_file);
q_beta_single = quaternion(transpose(beta_single_data));
% Estiamte beta_single ODF
beta_single_ori = orientation(q_beta_single, csCubic);
beta_single_odf = calcDensity(beta_single_ori,'kernel',deLaValleePoussinKernel,'halfwidth',5*degree);

%% Read beta_mixedphase data
beta_mixed_file = fopen(beta_mixed_oris);
beta_mixed_data = textscan(beta_mixed_file, '%f%f%f%f', 'HeaderLines', 1, 'CollectOutput', 1);
beta_mixed_data = beta_mixed_data{:};
beta_mixed_file = fclose(beta_mixed_file);
q_beta_mixed = quaternion(transpose(beta_mixed_data));
% Estiamte alpha_single ODF
beta_mixed_ori = orientation(q_beta_mixed, csCubic);
beta_mixed_odf = calcDensity(beta_mixed_ori,'kernel',deLaValleePoussinKernel,'halfwidth',5*degree);

%% plot all alpha
figure(1)
newMtexFigure('layout',[2,3], 'figSize','huge')
subplot(2,3,1)
plotPDF(alpha_single_odf, Miller(0, 0, 0, 2, csHex),'contourf', 'minmax')
title("(a) (0002)", 'FontSize', 26)
subplot(2,3,2)
plotPDF(alpha_single_odf, Miller(1, 0, -1, 0, csHex),'contourf', 'minmax');
title("(c) (10-10)", 'FontSize', 26)
subplot(2,3,3), plotPDF(alpha_single_odf, Miller(1, 1, -2, 0, csHex),'contourf', 'minmax')
title("(e) (11-20)", 'FontSize', 26)
subplot(2,3,4), plotPDF(alpha_mixed_odf, Miller(0, 0, 0, 2, csHex),'contourf', 'minmax')
title("(b) (0002)", 'FontSize', 26)
subplot(2,3,5), plotPDF(alpha_mixed_odf, Miller(1, 0, -1, 0, csHex),'contourf', 'minmax')
title("(d) (10-10)", 'FontSize', 26)
subplot(2,3,6), plotPDF(alpha_mixed_odf, Miller(1, 1, -2, 0, csHex),'contourf', 'minmax')
title("(f) (11-20)", 'FontSize', 26)
f = gcm;
for i=1:length(f.children)
    ax = getappdata(f.children(i),'sphericalPlot');
    ax.TL.FontSize=22; % TL= top left
    ax.BL.FontSize=22; % BL= bottom left
end
annotate(xvector,'label',{'RD'},'BackgroundColor','w','MarkerSize',15,'FontSize', 26)
annotate(yvector,'label',{'TD'},'BackgroundColor','w','MarkerSize',15,'FontSize', 26)
annotate(zvector,'label',{'ND'},'BackgroundColor','w','MarkerSize',15,'FontSize', 26)
setColorRange('equal')
mtexColorbar('title','mrd','FontSize', 26, 'location','eastoutside')
%suptitle('Alpha Single/Mixed phase Texture')
hold on
exportgraphics(gcf,'alpha_all.tif','Resolution',300)

%% plot all beta
figure(2)
newMtexFigure('layout',[2,3], 'figSize','huge')
subplot(2,3,1)
plotPDF(beta_single_odf, Miller(0, 0, 1, csCubic),'contourf', 'minmax')
title("(a) (001)", 'FontSize', 26)
subplot(2,3,2)
plotPDF(beta_single_odf, Miller(1, 1, 0, csCubic),'contourf', 'minmax')
title("(c) (110)", 'FontSize', 26)
subplot(2,3,3)
plotPDF(beta_single_odf, Miller(1, 1, 1, csCubic),'contourf', 'minmax')
title("(e) (111)", 'FontSize', 26)
subplot(2,3,4)
plotPDF(beta_mixed_odf, Miller(0, 0, 1, csCubic),'contourf', 'minmax')
title("(b) (001)", 'FontSize', 26)
subplot(2,3,5)
plotPDF(beta_mixed_odf, Miller(1, 1, 0, csCubic),'contourf', 'minmax')
title("(d) (110)", 'FontSize', 26)
subplot(2,3,6)
plotPDF(beta_mixed_odf, Miller(1, 1, 1, csCubic),'contourf', 'minmax')
title("(e) (111)", 'FontSize', 26)
f = gcm;
for i=1:length(f.children)
    ax = getappdata(f.children(i),'sphericalPlot');
    ax.TL.FontSize=22; % TL= top left
    ax.BL.FontSize=22; % TL= bottom left
end
annotate(xvector,'label',{'RD'},'BackgroundColor','w','MarkerSize',15,'FontSize', 26)
annotate(yvector,'label',{'TD'},'BackgroundColor','w','MarkerSize',15,'FontSize', 26)
annotate(zvector,'label',{'ND'},'BackgroundColor','w','MarkerSize',15,'FontSize', 26)
setColorRange('equal')
mtexColorbar('title','mrd','FontSize', 26, 'location','eastoutside')
%suptitle('Beta Single/Mixed phase Texture')
exportgraphics(gcf,'beta_all.tif','Resolution',300)

%% plot beta odfs
figure(3)
newMtexFigure( 'figSize','huge')
%subplot(1,2,1)
beta_single_odf.SS = specimenSymmetry('orthorhombic');
betasingle_pf = plotSection(beta_single_odf,'contourf','phi2',45*degree);
setColorRange('equal')
saveas(gcf, '../DAMASK_results/beta_single_odf.png')
%subplot(1,2,2)
figure(4)
beta_mixed_odf.SS = specimenSymmetry('orthorhombic');
plotSection(beta_mixed_odf,'contourf','phi2',45*degree);
setColorRange('equal')
mtexColorbar('title','mrd','FontSize', 36, 'location','eastoutside')
exportgraphics(gcf,'beta_odfs.png','Resolution',300)

disp('======== ALL DONE ========')