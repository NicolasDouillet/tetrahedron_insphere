%% tetrahedron_insphere
%
% Function to compute and display the insphere of a given tetrahedron.
%
% Author & support : nicolas.douillet (at) free.fr, 2022.
%
%
%% Syntax
%
% tetrahedron_insphere(A, B, C, D);
%
% tetrahedron_insphere(A, B, C, D, option_display);
%
% [I, r, rc] = tetrahedron_insphere(A, B, C, D, option_display);
%
%
%% Description
%
% tetrahedron_insphere(A, B, C, D) computes and displays the insphere of ABC tetrahedron.
%
% tetrahedron_insphere(A, B, C, D, option_display) displays ABCD tetrahedron and with its insphere when option_display is set either to
% logical true or real numeric 1, and doesn't when it is set to logical false or real numeric 0.
%
% [I, r, rc] = tetrahedron_insphere(A, B, C, D, option_display) stores the results in [I, r, rc] vector.
%
%
%% See also
%
% | <https://fr.mathworks.com/matlabcentral/fileexchange/119778-triangle-incircle-3d-2d?s_tid=srchtitle triangle incircle>
% | <https://fr.mathworks.com/matlabcentral/fileexchange/65574-tetrahedron-circumsphere?s_tid=srchtitle tetrahedron circumsphere> |
%
%
%% Input arguments
%
%        [Ax]
% - A = [Ay] : real column vector double. numel(A) = 3. One of the four ABCD vertices.
%        [Az]
%
%        [Bx]
% - B = [By] real column vector double. numel(A) = 3. One of the four ABCD vertices.
%        [Bz]
%
%        [Cx]
% - C = [Cy] real column vector double. numel(A) = 3. One of the four ABCD vertices.
%        [Cz]
%
%        [Dx]
% - D = [Dy] real column vector double. numel(A) = 3. One of the four ABCD vertices.
%        [Dz]
%
% - option_display : logical *true(1) / false(0), to enable/disable the display mode.
%
%
%% Output arguments
%
%        [Ix]
% - I = [Iy] : real column vector double. numel(I) = 3. The insphere centre.
%        [Iz]
%
% - r : real scalar double. The insphere radius.
%
% - rc : logical *true(1) / false(0). The return code. rc is true when the
%        outputs are valid and false when they are invalid (degenerated cases).
%
%
%% Example #1
% Random tetrahedron from the 3D space
V = 2*(rand(3,4)-0.5);
tetrahedron_insphere(V(:,1),V(:,2),V(:,3),V(:,4));

%% Example #2
% Regular tetrahedron in the unit ball
V1 = [0 0 1]';
V2 = [2*sqrt(2)/3 0 -1/3]';
V3 = [-sqrt(2)/3 sqrt(6)/3 -1/3]';
V4 = [-sqrt(2)/3 -sqrt(6)/3 -1/3]';
[I,radius,rc] = tetrahedron_insphere(V1,V2,V3,V4,true);
rc % expected : rc = 1

%% Example #3
% Flat / degenerated tetrahedron 
V1 = [0 0 0]';
V2 = [1 0 0]';
V3 = [0 1 0]';
V4 = [1 1 0]';
[I,radius,rc] = tetrahedron_insphere(V1,V2,V3,V4,true);
rc % expected : rc = 0