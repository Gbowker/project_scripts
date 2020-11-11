function [odfsection] = odfsection(odf, slice)
    
    odf.SS=specimenSymmetry('orthorhombic');
    plotSection(odf,'contourf','phi2',slice*degree);
    mtexColorbar('location','southoutside','FontSize', 20)
    
    disp('Drawn odfsection')
    
end