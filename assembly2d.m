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

function [K,C,M, KS, MS, R] = assembly2d(P, BulkElements, SurfElements)

Nbulk = length(P);
boundarynodes = unique(SurfElements(:));
Nsurf = length(boundarynodes);

% find first square element in mesh (they are all equal)
for i=1:length(BulkElements)
    if BulkElements(i).is_square
        Square = getLocalMatrices(copyElement2d(BulkElements(i)));
        MSq = Square.M(:);
        CSq = Square.C(:);
        KSq = Square.K(:);
        break
    end
end


% (START) MATRIX ASSEMBLY IN THE BULK
% Precompute the total number of entries needed for the bulk elements
nBulkElements = length(BulkElements);

% Initialize cell arrays to store the local results for each element
ii_cell = cell(nBulkElements, 1);
jj_cell = cell(nBulkElements, 1);
vm_cell = cell(nBulkElements, 1);
vc_cell = cell(nBulkElements, 1);
vk_cell = cell(nBulkElements, 1);

% Parallel loop for matrix assembly in the bulk
parfor i = 1:nBulkElements  % Each iteration processes one bulk element
    Element = BulkElements(i);
    eind = Element.Pind;  % Node indices for this element
    oind = ones(Element.NVert, 1);  % Vector of ones for Kronecker product
    
    % Local assembly of ii, jj, vm, vc, vk for this element
    ii_local = repmat(eind, length(oind), 1);
    jj_local = repelem(eind, length(oind));
    
    if Element.is_square
        vm_local = MSq;  % Use the precomputed square element mass matrix
        vc_local = CSq;  % Use the precomputed square element consistency matrix
        vk_local = KSq;  % Use the precomputed square element stiffness matrix
    else
        Element = getLocalMatrices(Element);  % Compute local matrices for this element
        vm_local = Element.M(:);  % Local mass matrix
        vc_local = Element.C(:);  % Local consistency matrix
        vk_local = Element.K(:);  % Local stiffness matrix
    end
    
    % Store the local results in the corresponding cell arrays
    ii_cell{i} = ii_local;
    jj_cell{i} = jj_local;
    vm_cell{i} = vm_local;
    vc_cell{i} = vc_local;
    vk_cell{i} = vk_local;
end

% After the parfor loop, concatenate the results from all elements
ii = vertcat(ii_cell{:});
jj = vertcat(jj_cell{:});
vm = vertcat(vm_cell{:});
vc = vertcat(vc_cell{:});
vk = vertcat(vk_cell{:});
% (END) MATRIX ASSEMBLY IN THE BULK


% (START) MATRIX ASSEMBLY ON THE SURFACE
% Precompute the total number of surface elements
nSurfElements = length(SurfElements);

% Precompute the lengths outside the loop
element_lengths = arrayfun(@(i) norm(P(SurfElements(i, 1), :) - P(SurfElements(i, 2), :)), 1:length(SurfElements));

% Initialize cell arrays to store results for each iteration
sii_cell = cell(nSurfElements, 1);
sjj_cell = cell(nSurfElements, 1);
svm_cell = cell(nSurfElements, 1);
svk_cell = cell(nSurfElements, 1);

% Parallel loop for assembling the matrix on the surface
parfor i = 1:nSurfElements
    % Extract the surface element indices
    eind = SurfElements(i, :);
    oind = [1; 1];
    
    % Calculate the local contributions for sii and sjj
    sii_local = repmat(eind, length(oind), 1);
    sjj_local = repelem(eind, length(oind));
    
    % Calculate the element length
    element_length = element_lengths(i);
    
    % Calculate the local contributions for svm and svk
    svm_local = [2 1 1 2]' / 6 * element_length;
    svk_local = [1 -1 -1 1]' / element_length;
    
    % Store local results in cell arrays
    sii_cell{i} = sii_local;
    sjj_cell{i} = sjj_local;
    svm_cell{i} = svm_local;
    svk_cell{i} = svk_local;
end

% Concatenate the results from the cell arrays into final vectors
sii = vertcat(sii_cell{:});
sjj = vertcat(sjj_cell{:});
svm = vertcat(svm_cell{:});
svk = vertcat(svk_cell{:});
% (END) MATRIX ASSEMBLY ON THE SURFACE

MS = sparse(sii,sjj, svm);
KS = sparse(sii,sjj, svk);
MS = MS(boundarynodes,boundarynodes);
KS = KS(boundarynodes,boundarynodes);
            
R = spalloc(Nbulk, Nsurf, Nsurf);
R(boundarynodes,:) = speye(Nsurf);

M = sparse(ii,jj, vm);
C = sparse(ii,jj, vc);
K = sparse(ii,jj, vk);

end