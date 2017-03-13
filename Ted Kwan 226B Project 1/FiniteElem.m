%% Finite Element Method
%
%   Written by Ted Kwan for Math 226B
%
%   Returns u(x,t) to approximate the solution
%   to the Laplace equation -\Delta u=f
%
%%% Inputs
%
% * mesht - String for mesh type (square or circle).
% * bdryt - String for boundary type (Dirichlet, Neumann or Mixed).
% * f - function handle for the right hand side.
% * ref - number of refinements to be used.
% * varargin options:
% 1 function handle for g(x).
% 2 function handle for g_n (normal derivative).
% * For mesht='circle'
% 3 x coordinate.
% 4 y coordinate.
% 5 radius r.
% 6 Mesh size h.
% * For mesht='square'
% 3 vector with the four boundaries of the rectangle.
% 4 Mesh size h.
%
function [u,node,elem,A]= FiniteElem(mesht,bdryt,f,ref,varargin)
    d=2; % Sets Dimension to two.
    %% Mesh Generation
    %
    %   Chooses a type of mesh generation based on
    %   the input parameter mesht.
    %
    if(strcmpi(mesht,'circle'))
        if(nargin<6) % Safety for access of elements.
            error('Not enough arguments');
        else
            %%% Create circle mesh
            %
            x=varargin{3};y=varargin{4};
            r=varargin{5};h=varargin{6};
            [node,elem] = circlemesh(x,y,r,h);
        end
    else
        if(nargin<4) % Safety for access of elements.
            error('Not enough input arguments');
        else
            %%% Create square mesh.
            %
            squ=varargin{3}; h=varargin{4};
            [node,elem] = squaremesh(squ,h);
        end
    end
    %%% Mesh Refinement.
    %
    while(ref>0)
        [node,elem] = uniformrefine(node,elem);
        ref=ref-1;
    end
    %% Assemble Stiffness Matrix.
    %
    %   Quick method to generate sparse matrix
    %   A, the stiffness matrix.
    %
    N=size(node,1); NT=size(elem,1);
    ii=zeros(9*NT,1); jj=zeros(9*NT,1); sA=zeros(9*NT,1);
    ve(:,:,3)=node(elem(:,2),:)-node(elem(:,1),:);
    ve(:,:,1)=node(elem(:,3),:)-node(elem(:,2),:);
    ve(:,:,2)=node(elem(:,1),:)-node(elem(:,3),:);
    area=0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));
    index=0;
    for i=1:3
       for j=1:3
           ii(index+1:index+NT)=elem(:,i);
           jj(index+1:index+NT)=elem(:,j);
           sA(index+1:index+NT)=dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
           index=index+NT;
       end
    end
    A=sparse(ii,jj,sA,N,N);
    %%% End Stiffness
    %% Calculate RHS
    %
    mid1 = (node(elem(:,2),:)+node(elem(:,3),:))/2;
    mid2 = (node(elem(:,3),:)+node(elem(:,1),:))/2;
    mid3 = (node(elem(:,1),:)+node(elem(:,2),:))/2;
    bt1 = area.*(f(mid2(:,1),mid2(:,2))+f(mid3(:,1),mid3(:,2)))/6;
    bt2 = area.*(f(mid3(:,1),mid3(:,2))+f(mid1(:,1),mid1(:,2)))/6;
    bt3 = area.*(f(mid1(:,1),mid1(:,2))+f(mid2(:,1),mid2(:,2)))/6;
    b = accumarray(elem(:),[bt1;bt2;bt3],[N,1]);
    %%% End RHS calculation.
    
    %% Boundary Condition.
    %
    if(strcmpi(bdryt,'Dirichlet'))
        %%% Dirichlet Boundary Conditions.
        %
        g_D=varargin{1};
        [bdNode,bdEdge,isBdNode]=findboundary(elem);
        freeNode = find(~isBdNode);
        u = zeros(N,1);
        u(bdNode) = g_D(node(bdNode,1),node(bdNode,2));
        b = b - A*u;
        u(freeNode)=A(freeNode,freeNode)\b(freeNode);
    elseif(strcmpi(bdryt,'Neumann'))
        %%% Neumann Boundary Conditions
        %
        g_N=varargin{2};
        [bdNode,bdEdge,isBdNode,isBdElem]=findboundary(elem);
        if(strcmpi(mesht,'circle'))
            bdFlag=SetBdFlagCir(isBdElem,elem,isBdNode,NT,d);
        else
            bdFlag=SetBdFlagSq(isBdElem,elem,isBdNode,bdEdge,NT,d);
        end
        u = zeros(N,1);
        totalEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
        Neumann = totalEdge(bdFlag(:) == 2,:);
        Nve = node(Neumann(:,1),:) - node(Neumann(:,2),:);
        edgeLength = sqrt(sum(Nve.^2,2));
        mid = (node(Neumann(:,1),:) + node(Neumann(:,2),:))/2;
        b = b + accumarray([Neumann(:),ones(2*size(Neumann,1),1)], ...
		repmat(edgeLength.*g_N(mid(:,1),mid(:,2))/2,2,1),[N,1]);
        b=b-mean(b)*ones(N,1); % Compatibility Condition.
        %u=bicgstabl(A,b,1e-8,500); %Uncomment to use 
        %u=gmres(A,b,3,1e-6,300);
        u=A\b;
        %%% Fixup for the additive constant.
        %
        %   This ensures that the constant is 0 as it
        %   is supposed to be, to compare it to the chosen
        %   real solution.
        ma=max(u);mi=min(u);
        if(ma>0)
            u=u-((ma)*ones(N,1));
        else
            u=u+(abs(ma)*ones(N,1));
        end
    else
        %%% Mixed boundary conditions - Not implemented
        %
        %   This method is not implemented, because the variable
        %   bdFlag is not passed to the function.
        %
        error('Mixed Boundary Conditions Are not yet implemented.');
        totalEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
        Dirichlet = totalEdge(bdFlag(:) == 1,:);
        Neumann = totalEdge(bdFlag(:) == 2,:);
        %---------- Dirichlet boundary conditions ----------
        isBdNode = false(N,1);
        isBdNode(Dirichlet) = true;
        bdNode = find(isBdNode);
        freeNode = find(~isBdNode);
        u = zeros(N,1);
        u(bdNode) = g_D(node(bdNode,1),node(bdNode,2));
        b = b - A*u;
        %---------- Neumann boundary conditions ----------
        if (~isempty(Neumann))
        Nve = node(Neumann(:,1),:) - node(Neumann(:,2),:);
        edgeLength = sqrt(sum(Nve.^2,2));
        mid = (node(Neumann(:,1),:) + node(Neumann(:,2),:))/2;
        b = b + accumarray([Neumann(:),ones(2*size(Neumann,1),1)], ...
        repmat(edgeLength.*g_N(mid(:,1),mid(:,2))/2,2,1),[N,1]);
        end
    end
end


