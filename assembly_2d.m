%
% DESCRIPTION - Assembles bulk and surface VEM meshes given a polygonal
% mesh
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS:
%
%  - P: array of nodes
%  - BulkElements: all elements of the bulk
%  - SurfaceElements: all elements of the surface
%
% OUTPUTS:
%
% - K,C,M: stiffness, consistency and mass matrices in the bulk
% - KS,MS: stiffness and mass matrices on the surface
% - R: reduction matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K,C,M, KS, MS, R] = assembly_2d(P, BulkElements, SurfElements)

Nbulk = length(P);
boundarynodes = unique(SurfElements(:));
Nsurf = length(boundarynodes);

K = spalloc(Nbulk,Nbulk,9*Nbulk); % Stiffness matrix in the 2D bulk
C = spalloc(Nbulk,Nbulk,9*Nbulk); % Mass matrix in the 2D bulk
M = spalloc(Nbulk,Nbulk,9*Nbulk); % Mass matrix in the 2D bulk

KS = spalloc(Nbulk, Nbulk, 3*Nsurf); % Stiffness matrix on the 1D surface
MS = spalloc(Nbulk, Nbulk, 3*Nsurf); % Mass matrix on the 1D surface

% find first square element in mesh (they are all equal)
for i=1:length(BulkElements)
    if BulkElements(i).is_square
        Square = copyElement2d(BulkElements(i));
        MSq = getLocalMatrices(Square).M;
        CSq = getLocalMatrices(Square).C;
        KSq = getLocalMatrices(Square).K;
        break
    end
end

% MATRIX ASSEMBLY IN THE BULK
for i=1:length(BulkElements) % For each bulk element
    Element = copyElement2d(BulkElements(i));
    eind = Element.Pind;
    if Element.is_square
        M(eind, eind) = M(eind, eind) + MSq; %#ok
        C(eind, eind) = C(eind, eind) + CSq; %#ok
        K(eind, eind) = K(eind, eind) + KSq; %#ok
    else
        M(eind, eind) = M(eind, eind) + getLocalMatrices(Element).M; %#ok
        C(eind, eind) = C(eind, eind) + getLocalMatrices(Element).C; %#ok
        K(eind, eind) = K(eind, eind) + getLocalMatrices(Element).K; %#ok
    end
end

% MATRIX ASSEMBLY ON THE SURFACE
for i=1:length(SurfElements)
    eind = SurfElements(i,:);
    element_length = norm(P(eind(1),:) - P(eind(2),:));
    MS(eind, eind) = MS(eind, eind) + [2 1; 1 2]/6*element_length; %#ok
    MS(eind, eind) = MS(eind, eind) + [1 -1; -1 1]/element_length; %#ok
end

MS = MS(boundarynodes,boundarynodes);
KS = KS(boundarynodes,boundarynodes);
            
R = spalloc(Nbulk, Nsurf, Nsurf);
R(boundarynodes,:) = speye(Nsurf);

end