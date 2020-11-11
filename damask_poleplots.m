%% Define MTEX preferences...
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outofPlane');
setMTEXpref('figSize','large');
setMTEXpref('FontSize', 20);

setMTEXpref('defaultColorMap', 'Viridis');

%% === load ti_phases: === %

% --- User defined --- %
phase = 'alpha';                                                            % change to 'alpha' or 'beta'
miller_indicies = [0, 0, 0, 2];                                             % change to m or mb
% -------------------- %

[cs, h, k, l, w, m_str] = ti_phases(phase, miller_indicies);                % function returns phase parameters

%% Define path to quaternion.txt files...

% --- User defined --- %
path_to_oris = '../DAMASK_results/mixed_1024grains/alpha/';                 % define path to containing dir
ori_increment = {'ori_initial.txt', 'ori_final.txt'};                       % define files
% -------------------- %

%% for each increment file...

for index = 1:length(ori_increment)
    
    full_path = strcat(path_to_oris, ori_increment{index});                 % define the file
    
    % Read the quaternions from the file
    fid = fopen(full_path);
    data = textscan(fid, '%f%f%f%f', 'HeaderLines', 1, 'CollectOutput', 1);
    data = data{:};
    fid = fclose(fid);
    q = quaternion(transpose(data));
    
    
    % Estiamte an ODF from the orientations
    ori = orientation(q, cs);
    
    % calcdensity method
    odf = calcDensity(ori,'kernel',deLaValleePoussinKernel,'halfwidth',10*degree);
 
    % calcKernel method
    %psi = calcKernel(ori,'Method','KLCV')
    %odf = calcKernelODF(ori,'kernel',psi);
    
    % rotating results
    odf = rotate(odf,rotation.byAxisAngle(yvector,0*degree));
    
    % plot pole figures
    figure(index)
    poleplot(odf, phase, h, k, l, w, cs, m_str)
    exportgraphics(gcf, strcat(path_to_oris, './pdf_',m_str, '_', ori_increment{index},'.tiff'), 'Resolution', 200)
    
    % plot odf phi2 slices for beta
    if isequal(phase, 'beta')
        figure(index+2)
        odfsection(odf, 45)
        exportgraphics(gcf, strcat(path_to_oris, './odf_', ori_increment{index}, '.tiff'), 'Resolution', 200)
    end
end

disp('======== ALL DONE ========')