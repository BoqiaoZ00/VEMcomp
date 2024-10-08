function [P, h, BulkElements, SurfElements] = generate_mesh2d(fun, Q, Nx, tol)
%GENERATE_MESH_FLAT_DOMAIN Summary of this function goes here
%   Detailed explanation goes here

hx_requested = (Q(1,2)-Q(1,1))/(Nx-1); % Requested discretisation step along x
hy_requested = (Q(2,2)-Q(2,1))/(Nx-1); % Requested discretisation step along y

h_mono = min([hx_requested,hy_requested]); % Actual discretization step along each dimension
h = h_mono*sqrt(2); % Meshsize

Nx = ceil((Q(1,2)-Q(1,1))/h_mono)+1; % Corrected number of discretization nodes along x
Ny = ceil((Q(2,2)-Q(2,1))/h_mono)+1; % Corrected number of discretization nodes along y
Nrect = Nx*Ny; % Number of nodes of bounding box

Q(1,:) = Q(1,:) + (h_mono*(Nx-1) - (Q(1,2)-Q(1,1)))/2*[-1,1]; % Corrected range of bounding box along x
Q(2,:) = Q(2,:) + (h_mono*(Ny-1) - (Q(2,2)-Q(2,1)))/2*[-1,1]; % Corrected range of bounding box along y

% GENERATE ONE SQUARE ELEMENT
PS = [0 0 0; 0 1 0; 1 1 0; 1 0 0]*h_mono;
ES = element2d(PS, true);

% (CODE MODIFIED) CREATING GRIDPOINTS OF BOUNDING BOX
[X, Y] = meshgrid(linspace(Q(1,1), Q(1,2), Nx), linspace(Q(2,1), Q(2,2), Nx));
P = [X(:), Y(:), zeros(Nx^2, 1)];

% (CODE MODIFIED) GENERATE ELEMENTS
max_num_of_elements = (Nx-1)*(Ny-1);
count_square_ele = 0;
count_nonsquare_ele = 0;
count_se = 0;
dummy_element2d = element2d();
SquareElements(max_num_of_elements, 1) = dummy_element2d;
NonSquareElements(max_num_of_elements, 1) = dummy_element2d;
SE = repmat({zeros(2, 3)}, max_num_of_elements, 1); 
accepted_node = false(Nrect,1);
for i=0:Nx-2 % For each element of the bounding box
    for j=0:Ny-2
        indexes = [Ny*i+j+1
                   Ny*i+j+2
                   Ny*(i+1)+j+2
                   Ny*(i+1)+j+1];                  
        NewSquareElement = shiftElement(ES, P(indexes(1),:));
        if is_outside(NewSquareElement, fun)
            continue
        end
        if is_inside(NewSquareElement, fun)
            % Store square elements that are inside domain
            count_square_ele = count_square_ele + 1;
            accepted_node(indexes) = true(4,1);
            setPind(NewSquareElement, indexes);
            SquareElements(count_square_ele) = NewSquareElement; 
            continue
        end
        % Store non-square elements obtained by cutting square elements
        % with boundary. Such non-square elements are not endowed with node
        % indexes yet.
        count_nonsquare_ele = count_nonsquare_ele + 1;
        count_se = count_se + 1;
        [NewElement, LocalSurfaceElements] = cut(NewSquareElement, fun, tol);
        NonSquareElements(count_nonsquare_ele) = NewElement; 
        SE(count_se) = LocalSurfaceElements;
    end
end
SquareElements = SquareElements(1:count_square_ele); 
NonSquareElements = NonSquareElements(1:count_nonsquare_ele); 
SE = SE(1:count_se); 

% AFTER ELIMINATING SQUARE ELEMENTS THAT ARE OUTSIDE DOMAIN, RE-DETERMINE
% INDEXES OF NODES USED BY SQUARE ELEMENTS
P = P(accepted_node,:);
acceptedindexes = zeros(Nrect,1);
acceptedindexes(accepted_node,1) = linspace(1,length(P),length(P))';
for i=1:length(SquareElements)
    setPind(SquareElements(i), acceptedindexes(SquareElements(i).Pind));
end

% DETERMINE SET OF NON-REPEATED NODES UP TO SMALL TOLERANCE
for i=1:length(NonSquareElements)
   P = [P; NonSquareElements(i).P]; %#ok 
end
P = uniquetol(P,tol,'ByRows',true);

% FIX ELEMENTS BY ASSIGNING NODE INDEXES AND ELIMINATING DUPLICATE NODES UP
% TO SMALL TOLERANCE

% QUESTION: if two SquareElements are with the tol, then they refer to the
% same element in P. But in this case, does SquareElements list contain
% duplicate elements?
for i=1:length(SquareElements)
   [~, ind] = ismembertol(SquareElements(i).P,P,tol,'ByRows',true);
   setPind(SquareElements(i), ind);
   setP(SquareElements(i), P(ind,:));
end

for i=1:length(NonSquareElements)
   [~, ind] = ismembertol(NonSquareElements(i).P,P,tol,'ByRows',true);
   setPind(NonSquareElements(i), ind);
   setP(NonSquareElements(i), P(ind,:));
end

BulkElements = [SquareElements; NonSquareElements];

SurfElements = zeros(length(SE),2);
for i=1:length(SE)
    [~, ind] = ismembertol(SE{i},P,tol,'ByRows',true);
    SurfElements(i,:) = ind';
end

end

% DETERMINE IF ELEMENT IS OUTSIDE DOMAIN (POSSIBLY TOUCHING BOUNDARY)
function outside = is_outside(Element, fun)
    values = fun(Element.P);
    outside = all(values >= 0);
end

% DETERMINE IF ELEMENT IS INSIDE DOMAIN (POSSIBLY TOUCHING BOUNDARY)
function inside = is_inside(Element, fun)
    values = fun(Element.P);
    inside = all(values <= 0);
end

% CUTS GIVEN ELEMENT BY BOUNDARY OF DOMAIN
function [CutElement, LocalSurfaceElements] = cut(Element, fun, tol)
    P = Element.P;
    inside_or_boundary = fun(P) <= 0;
    Pnew = [];
    for i=1:length(P)
       j = rem(i, length(P)) + 1;
       if inside_or_boundary(i)
           Pnew = [Pnew; P(i,:)]; %#ok
       end
       if sum(inside_or_boundary([i j])) == 1
           Pnew = [Pnew; intersectEdge(P(i,:), P(j,:), fun)]; %#ok
       end
    end
    [~, ind] = uniquetol(Pnew,tol,'Byrows',true);
    if length(ind) < 3
        CutElement = [];
        return
    end
    indnewsort = sort(ind);
    Pnewsort = Pnew(indnewsort,:);
    CutElement = element2d(Pnewsort, false, false, [], mean(Pnewsort,1));
    indboundary = find(abs(fun(Pnewsort)) < tol);
    LocalSurfaceElements = cell(0,1);
    for i=1:length(indnewsort)-1
        if all(ismember(indnewsort([i i+1]), indboundary))
            LocalSurfaceElements = [LocalSurfaceElements; {Pnew(indnewsort([i i+1]),:)}]; %#ok
        end
    end
    if all(ismember(indnewsort([end 1]), indboundary))
            LocalSurfaceElements = [LocalSurfaceElements; {Pnew(indnewsort([end 1]),:)}];
    end
end

% COMPUTES INTERSECTION POINT BETWEEN AN EDGE AND THE BOUNDARY OF THE
% DOMAIN
function newP = intersectEdge(P1, P2, fun)
    fun_restricted = @(alpha) fun(alpha*P1 + (1-alpha)*P2);
    alpha = fzero(fun_restricted, [0 1]);
    newP = alpha*P1 + (1-alpha)*P2;
end