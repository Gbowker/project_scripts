function [cs, h, k, l, w, m_str] = ti_phases(phase, miller_indicies)
    
    h = miller_indicies(1);
    k = miller_indicies(2);
    l = miller_indicies(3);
    m_str = strcat(string(h), string(k), string(l));
    
    if isequal(phase,'alpha')   % alpha phase %
        cs = crystalSymmetry('6/mmm', [3 3 4.7], 'X||b', 'Y||a*', 'Z||c*', 'mineral', 'Ti-Hex', 'color', [0.53 0.81 0.98]);
        %cs = crystalSymmetry('6/mmm', [3 3 4.7], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Ti-Hex', 'color', [0.53 0.81 0.98]);
        w = miller_indicies(4);
        m_str = strcat(m_str, string(w));
    else                        % beta phase %
        cs = crystalSymmetry('m-3m', [3.2 3.2 3.2], 'mineral', 'Titanium cubic', 'color', [0.56 0.74 0.56]);
        w = 0;
    end
    disp('got phase paramters...')
end

% Notes: %
% crystalSymmetry('laue group', [lattice parameters(unit cell size)], axes convention, 'mineral', 'ebsddatacolor', [R G B]);