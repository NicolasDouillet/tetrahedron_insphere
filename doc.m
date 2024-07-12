%% tetrahedron_insphere
%
% Function to compute and display the insphere of a given tetrahedron.
%
% Author : nicolas.douillet (at) free.fr, 2022-2024.
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
% tetrahedron_insphere(A, B, C, D) computes and displays the insphere of ABCD tetrahedron.
%
% tetrahedron_insphere(A, B, C, D, option_display) displays ABCD tetrahedron and with its insphere when
% option_display is set either to logical true or real numeric 1, and doesn't when it is set to logical
% false or real numeric 0.
%
% [I, r, rc] = tetrahedron_insphere(A, B, C, D, option_display) stores the results in [I, r, rc] vector.
%
%
%% See also
%
% | <https://fr.mathworks.com/matlabcentral/fileexchange/119778-triangle-incircle-3d-2d?s_tid=srchtitle triangle incircle> |
%   <https://fr.mathworks.com/help/matlab/ref/triangulation.incenter.html incenter> |
%   <https://fr.mathworks.com/matlabcentral/fileexchange/65574-tetrahedron-circumsphere?s_tid=srchtitle tetrahedron circumsphere> |
%
%
%% Input arguments
%
%        [Ax]
% - A = [Ay] : real column vector double. numel(A) = 3. One of the four ABCD vertices.
%        [Az]
%
% - B, C, D : same type and description as A, here above.
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
A = [0 0 1]';
B = [2*sqrt(2)/3 0 -1/3]';
C = [-sqrt(2)/3 sqrt(6)/3 -1/3]';
D = [-sqrt(2)/3 -sqrt(6)/3 -1/3]';
[I,r] = tetrahedron_insphere(A,B,C,D,true) % expected : I = [0 0 0]; r = 1/3;

%% Example #3
% Flat / degenerated tetrahedron 
A = [0 0 0]';
B = [1 0 0]';
C = [0 1 0]';
D = [1 1 0]';
[I,r,rc] = tetrahedron_insphere(A,B,C,D,true);
rc % expected : rc = 0