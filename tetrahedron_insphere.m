function [I, r, rc] = tetrahedron_insphere(A, B, C, D, option_display)
%% tetrahedron_insphere : function to compute and display the insphere of a given tetrahedron.
%
% Author & support : nicolas.douillet (at) free.fr, 2022.
%
%
% Syntax
%
% tetrahedron_insphere(A, B, C, D);
% tetrahedron_insphere(A, B, C, D, option_display);
% [I, r, rc] = tetrahedron_insphere(A, B, C, D, option_display);
%
%
% Description
%
% tetrahedron_insphere(A, B, C, D) computes and displays the insphere of ABC tetrahedron.
% tetrahedron_insphere(A, B, C, D, option_display) displays ABCD tetrahedron and with its insphere when option_display is set either to
% logical true or real numeric 1, and doesn't when it is set to logical false or real numeric 0.
% [I, r, rc] = tetrahedron_insphere(A, B, C, D, option_display) stores the results in [I, r, rc] vector.
%
%
% See also SPHERE
%
%
% Input arguments
%
%       [Ax]
% - A = [Ay] : real column vector double. numel(A) = 3. One of the four ABCD vertices.
%       [Az]
%
%       [Bx]
% - B = [By] real column vector double. numel(A) = 3. One of the four ABCD vertices.
%       [Bz]
%
%       [Cx]
% - C = [Cy] real column vector double. numel(A) = 3. One of the four ABCD vertices.
%       [Cz]
%
%       [Dx]
% - D = [Dy] real column vector double. numel(A) = 3. One of the four ABCD vertices.
%       [Dz]
%
% - option_display : logical *true(1) / false(0), to enable/disable the display mode.
%
%
% Output arguments
%
%       [Ix]
% - I = [Iy] : real column vector double. numel(I) = 3. The insphere centre.
%       [Iz]
%
% - r : real scalar double. The insphere radius.
%
% - rc : logical *true(1) / false(0). The return code. rc is true when the
%        outputs are valid and false when they are invalid (degenerated cases).
%
%
% Example #1
% Random tetrahedron from the 3D space
% V = 2*(rand(3,4)-0.5);
% [I,radius] = tetrahedron_insphere(V(:,1),V(:,2),V(:,3),V(:,4),true);
%
%
% Example #2
% Regular tetrahedron in the unit ball
% V1 = [0 0 1]';
% V2 = [2*sqrt(2)/3 0 -1/3]';
% V3 = [-sqrt(2)/3 sqrt(6)/3 -1/3]';
% V4 = [-sqrt(2)/3 -sqrt(6)/3 -1/3]';
% [I,radius,rc] = tetrahedron_insphere(V1,V2,V3,V4,true);
% rc % expected : rc = 1
%
%
% Example #3
% Flat / degenerated tetrahedron 
% V1 = [0 0 0]';
% V2 = [1 0 0]';
% V3 = [0 1 0]';
% V4 = [1 1 0]';
% [I,radius,rc] = tetrahedron_insphere(V1,V2,V3,V4,true);
% rc % expected : rc = 0


%% Input parsing
assert(nargin > 3, 'Not enought input arguments. Three input points are required to define one tetrahedron.');
assert(nargin < 6, 'Too many input arguments.');

if nargin < 5    
    option_display = true;            
end

assert(isequal(size(A),size(B),size(C),size(D)),'All inputs points must have the same size.');
assert(isequal(ndims(A),ndims(B),ndims(C),ndims(D),2),'All inputs points must have the same number of dimensions (2).');
assert(isreal(A) && isreal(B) && isreal(C) && isreal(D),'All inputs points must contain real numbers only.');
assert(numel(A) == 3,'Input points must have exactly 3 elements.');


%% Body

% 0 Process ABCD degenerated cases (flat or aligned points)
rc = true;
dist_D_2_ABC = point_to_plane_distance(D',cross(B'-A',C'-A'),A');

if dist_D_2_ABC < 1e4*eps    
    warning('ABCD tetrahedron is flat or degenerated; ABCD insphere and its centre are irrelevant.');
    rc = false;
end

% I Compute ABCD tetrahedron edge vectors
AB = (B-A);
AC = (C-A);
AD = (D-A);
BC = (C-B);

% II.1 Project B and C on AD dihedron
[~,H_B_AD] = point_to_plane_distance(B',AD',A');
[~,H_C_AD] = point_to_plane_distance(C',AD',A');

% II.2 Project C and D on AB dihedron
[~,H_C_AB] = point_to_plane_distance(C',AB',A');
[~,H_D_AB] = point_to_plane_distance(D',AB',A');

% II.3 Project A and D on BC dihedron
[~,H_A_BC] = point_to_plane_distance(A',BC',B');
[~,H_D_BC] = point_to_plane_distance(D',BC',B');

% III Project these points on the plane orthogonal to the dihedron
% axis and going through the previous angle point. Post normalization step
% is essential for the following.
H_B_AD = (H_B_AD - A') / norm((H_B_AD - A')) + A';
H_C_AD = (H_C_AD - A') / norm((H_C_AD - A')) + A';
H_C_AB = (H_C_AB - A') / norm((H_C_AB - A')) + A';
H_D_AB = (H_D_AB - A') / norm((H_D_AB - A')) + A';
H_A_BC = (H_A_BC - B') / norm((H_A_BC - B')) + B';
H_D_BC = (H_D_BC - B') / norm((H_D_BC - B')) + B';

% IV Compute the average points, belonging to the dihedra bisector planes.
H_AD = 0.5 * (H_B_AD + H_C_AD)';
H_AB = 0.5 * (H_C_AB + H_D_AB)';
H_BC = 0.5 * (H_A_BC + H_D_BC)';

% V Compute the corresponding diedral planes (each one of them
% defined by three points, the two of the diedral axis, plus the one
% projected on its orthogonal plane.
n1 = cross((H_AD-A),AD);
n4 = cross((H_AB-A),AB);
n5 = cross((H_BC-B),BC);

% VI Compute the intersection of two dihedral bisector planes. Result is
% the line of the 3D space.
[I1,u1] = planes_intersection(n1,H_AD,n4,H_AB);

% VII Compute intersection between this line and one third dihedral
% bisector plane; this is ABCD tetrahedron insphere centre.
I = line_plane_intersection(u1,I1,n5,H_BC);

% VIII Compute the distance between this point and one of the four
% tetrahedron faces; this is ABCD tetrahedron insphere radius.
n = cross(AB,AC);
r = point_to_plane_distance(I',n',A');

% Sphere
[Sx,Sy,Sz] = sphere(60);


%% Display
if option_display
    
    V = cat(2,A,B,C,D);
    cmap_tetra = cat(3,ones(size(Sx)),zeros(size(Sy)),zeros(size(Sz)));
    vtx_triplets = combnk(1:4,3);
    
    figure;
    plot3(I(1,1),I(2,1),I(3,1),'ko','Linewidth',4), hold on;
    s = surf(r*Sx+I(1,1),r*Sy+I(2,1),r*Sz+I(3,1),cmap_tetra); shading interp;
    alpha(s,0.4);
    
    for k = 1:size(vtx_triplets,1)
        f = fill3(V(1,vtx_triplets(k,:)),V(2,vtx_triplets(k,:)),V(3,vtx_triplets(k,:)),'b','EdgeColor','b'); hold on;
        alpha(f,0.2);
    end
    
    axis equal, axis tight;
    ax = gca;
    ax.Clipping = 'off';
    [az,el] = view(3);
    camlight(az,el);
    
end


end % tetrahedron_insphere


%% point_to_plane_distance subfunction
function [d2H, H] = point_to_plane_distance(M, n, I)
%
% Author and support nicolas.douillet (at) free.fr, 2020.

nb_pts = size(M,1);
d_I = -(n(1)*I(1)+n(2)*I(2)+n(3)*I(3));
t_H = -(repmat(d_I,[nb_pts,1])+n(1)*M(:,1)+n(2)*M(:,2)+n(3)*M(:,3)) / sum(n.^2);

x_H = M(:,1) + t_H*n(1);
y_H = M(:,2) + t_H*n(2);
z_H = M(:,3) + t_H*n(3);

% Orthogonal projected point
H = cat(2,x_H,y_H,z_H);
d2H = sqrt(sum((M-H).^2,2));


end % point_to_plane_distance


%% planes_intersection subfunction
function [I, u, rc] = planes_intersection(n1, M1, n2, M2)
%
% Author & support : nicolas.douillet (at) free.fr, 2019-2022.

d1 = -dot(n1,M1); % -a1*x1 - b1*y1 - c1*z1
d2 = -dot(n2,M2); % -a2*x2 - b2*y2 - c2*z2

u = cross(n1,n2);

if norm(u) == 0 % (M1,n1) = (M2,n2) or (M1,n1) // (M2,n2)
   
    if (dot(n1,M2) + d1) == 0 && (dot(n2,M1) + d2) == 0 % (a1*M2(1) + b1*M2(2) + c1*M2(3) + d1) == 0
                        
        I = M1;
        u = M2 - M1;
        rc = 2;
        
    else                
        
        I = [];
        u = [];
        rc = 0;
        
    end
    
else 
          
     dir = find((abs(u) == max(abs(u))));     
     dir = dir(1);
     
     % => the line does exist in this direction, and then it can be set to t = 0.
     
     switch dir
         
         case 1 % setting : x = 0
             
             dx0y = (n1(3)*d2 - n2(3)*d1); % c1*d2 - c2*d1
             dx0z = (n2(2)*d1 - n1(2)*d2); % b2*d1 - b1*d2
             
             xI = 0;           
             yI = dx0y/u(1); 
             zI = dx0z/u(1);
             
         case 2 % setting : y = 0
             
             dxy0 = (n1(3)*d2 - n2(3)*d1); % c1*d2 - c2*d1
             dy0z = (n2(1)*d1 - n1(1)*d2); % a2*d1 - a1*d2
             
             xI = -dxy0/u(2);
             yI = 0;
             zI = -dy0z/u(2);
             
         case 3 % setting : z = 0
             
             dxz0 = (n1(2)*d2 - n2(2)*d1); % b1*d2 - b2*d1
             dyz0 = (n2(1)*d1 - n1(1)*d2); % a2*d1 - a1*d2
             
             xI = dxz0/u(3);
             yI = dyz0/u(3);
             zI = 0;                         
             
     end
     
     I = zeros(size(M1));
     I(1) = xI;
     I(2) = yI;
     I(3) = zI;
     
     rc = 1;
     
end


end % planes_intersection


%% line_plane_intersection subfunction
function [I,rc] = line_plane_intersection(u, N, n, M)
%
% Author & support : nicolas.douillet (at) free.fr, 2020 - 2022.

% Plane offset parameter
d = -dot(n,M);

% Specific cases treatment
if ~dot(n,u) % n & u perpendicular vectors
    if dot(n,N) + d == 0 % N in P => line belongs to the plane        
        I = N;
        rc = 2;
    else % line // to the plane        
        I = [];
        rc = 0;
    end
else
    
    % Parametric line parameter t
    t = - (d + dot(n,N)) / dot(n,u);
    
    % Intersection coordinates
    I = N + u*t;
    
    rc = 1;
    
end


end % line_plane_intersection