%
% DESCRIPTION - Assembles bulk and surface VEM meshes given a polyhedral
% mesh
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS:
%
%  - P: array of nodes
%  - Elements: all elements of the bulk
%  - EGamma: all elements of the surface
%
% OUTPUTS:
%
% - K,M,C: stiffness, mass, and consistency matrices in the bulk
% - KS,MS,CMS: stiffness, mass, and consistency matrices on the surface
% - R: reduction matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K, C, M, KS, CS, MS, R] = assembly3d(P, BulkElements, SurfElements)

Nbulk = length(P);
boundarynodes = unique(SurfElements(:));
Nsurf = length(boundarynodes);

if isempty(BulkElements) % Triangulated surface mesh
    K = []; M = []; C = []; R = speye(Nsurf);
    [KS,MS] = matrices(P,SurfElements);
    CS = MS;
    return
end

% find first cubic element in mesh (they are all equal)
for i=1:length(BulkElements)
    if BulkElements(i).is_cube
        Cube = getLocalMatrices(copyElement3d(BulkElements(i)));
        MC = Cube.M(:);
        KC = Cube.K(:);
        CC = Cube.C(:);
        % An element3dcube is not supposed to have boundary faces
        break
    end
end

% allocate vectors for sparse matrix representation
Nbulk_repeated = 0;
for i=1:length(BulkElements)
    Nbulk_repeated = Nbulk_repeated + BulkElements(i).NVert;
end
% Preallocate arrays with initial size (considering maximum possible size)
ii_local = cell(length(BulkElements), 1);
jj_local = cell(length(BulkElements), 1);
vk_local = cell(length(BulkElements), 1);
vm_local = cell(length(BulkElements), 1);
vc_local = cell(length(BulkElements), 1);

sii_local = cell(length(BulkElements), 1);
sjj_local = cell(length(BulkElements), 1);
svk_local = cell(length(BulkElements), 1);
svm_local = cell(length(BulkElements), 1);
svc_local = cell(length(BulkElements), 1);

% Parallel loop
parfor (i = 1:length(BulkElements), 4)  % For each bulk element
    E = BulkElements(i);
    eind = E.Pind;
    nVert = E.NVert;
    
    % Local storage for current element
    ii_local{i} = repmat(eind, nVert, 1);
    jj_local{i} = reshape(repmat(eind', nVert, 1), [], 1);

    if E.is_cube
        vm_local{i} = MC;
        vc_local{i} = CC;
        vk_local{i} = KC;
    else
        Element = getLocalMatrices(copyElement3d(E));
        vm_local{i} = Element.M(:);
        vc_local{i} = Element.C(:);
        vk_local{i} = Element.K(:);

        % Local storage for surface elements
        totalRows = sum(arrayfun(@(E) sum([E.Faces([E.Faces.is_boundary]).NVert].^2), Element));
        sii_temp = zeros(totalRows, 1);
        sjj_temp = zeros(totalRows, 1);
        svm_temp = zeros(totalRows, 1);
        svc_temp = zeros(totalRows, 1);
        svk_temp = zeros(totalRows, 1);

        sii_temp_count = 0;
        sjj_temp_count = 0;
        svm_temp_count = 0;
        svc_temp_count = 0;
        svk_temp_count = 0;
        

        for j = 1:Element.NFaces
            Face = copyElement2d(Element.Faces(j));
            if Face.is_boundary
                eind_boundary = Face.Pind;
                n = Face.NVert;  % Number of vertices

                % Use repmat and reshape to replicate indices instead of kron
                sii_temp(sii_temp_count+1:sii_temp_count+n*n) = repmat(eind_boundary, n, 1);
                sii_temp_count = sii_temp_count + n*n;
                sjj_temp(sjj_temp_count+1:sjj_temp_count+n*n) = reshape(repmat(eind_boundary', n, 1), [], 1);
                sjj_temp_count = sjj_temp_count + n*n;

                Face = getLocalMatrices(Face);
                svm_temp(svm_temp_count+1:svm_temp_count+n*n) = Face.M(:);
                svm_temp_count = svm_temp_count + n*n;
                svc_temp(svc_temp_count+1:svc_temp_count+n*n) = Face.C(:);
                svc_temp_count = svc_temp_count + n*n;
                svk_temp(svk_temp_count+1:svk_temp_count+n*n) = Face.K(:);
                svk_temp_count = svk_temp_count + n*n;
            end
        end
        sii_local{i} = sii_temp;
        sjj_local{i} = sjj_temp;
        svm_local{i} = svm_temp;
        svc_local{i} = svc_temp;
        svk_local{i} = svk_temp;
    end
end
% Combine local arrays into global arrays

ii = vertcat(ii_local{:});
jj = vertcat(jj_local{:});
vk = vertcat(vk_local{:});
vm = vertcat(vm_local{:});
vc = vertcat(vc_local{:});

% Combine surface element data
sii = vertcat(sii_local{:});
sjj = vertcat(sjj_local{:});
svk = vertcat(svk_local{:});
svm = vertcat(svm_local{:});
svc = vertcat(svc_local{:});

common_indices = [sii, sjj];
values = [svm, svc, svk];  % Combine all values
MS = sparse(common_indices(:, 1), common_indices(:, 2), values(:, 1));
CS = sparse(common_indices(:, 1), common_indices(:, 2), values(:, 2));
KS = sparse(common_indices(:, 1), common_indices(:, 2), values(:, 3));


MS = MS(boundarynodes,boundarynodes);
CS = CS(boundarynodes,boundarynodes);
KS = KS(boundarynodes,boundarynodes);

            
R = spalloc(Nbulk, Nsurf, Nsurf);
R(boundarynodes,:) = speye(Nsurf);

% Combine common parts in the matrix assembly (if suitable)
common_indices = [ii, jj];
values = [vm, vc, vk];  % Combine all values
M = sparse(common_indices(:, 1), common_indices(:, 2), values(:, 1));
C = sparse(common_indices(:, 1), common_indices(:, 2), values(:, 2));
K = sparse(common_indices(:, 1), common_indices(:, 2), values(:, 3));
        
end