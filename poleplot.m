function [poleplot] = poleplot(odf, phase, h, k, l, w, cs, m_str)
    
    if isequal(phase,'alpha')
        plotPDF(odf, Miller(h, k, l, w, cs), 'contourf', 'minmax')
        title(strcat('(', string(h), string(k), string(l), string(w), ')'), 'FontSize', 20)
    else
        plotPDF(odf, Miller(h, k, l, cs), 'contourf', 'minmax')
        title(strcat('(', string(h), string(k), string(l), ')'), 'FontSize', 26)
    end
    annotate(xvector,'label',{'RD'},'BackgroundColor','w','FontSize', 11)
    annotate(yvector,'label',{'TD'},'BackgroundColor','w','FontSize', 11)
    annotate(zvector,'label',{'ND'},'BackgroundColor','w','FontSize', 11)
    mtexColorbar('location','southoutside','FontSize', 20)
    
    f = gcm;
for i=1:length(f.children)
    ax = getappdata(f.children(i),'sphericalPlot');
    ax.TL.FontSize=20; % TL= top left
    ax.BL.FontSize=20; % BL= bottom left
end 
    disp('Drawn poleplot')
end