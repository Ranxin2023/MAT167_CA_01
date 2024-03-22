%% CA1
%Name: Randy Li
%Student ID:917196816
%Date:February 4th 2024
%load the file CA_01.mat
load CA_01;

%draw the signal x
%open a new figure called figure 1
figure(1);
stem(x); hold on; plot(x); grid;

%using discrete cosine transform to find orthonormal basis vectors 
fplot(@(x) cos(3*pi*(x-1)/7), [1,8], 'r--o');

%create a 8*1 matrix with all ones
A = ones(8, 1);
% loop 7 times to for the rest of the column vectors
for k = 1:7
    % loop to append column matrix by using the forluma in the instructions
    A = [A arrayfun(@(x) cos(k*pi*(x-1)/7), (1:8)')];

end

%display matrix A
fprintf("The matrix A is:\n")
disp(A)

%display ATA
fprintf('The result of ATA is\n');
ATA = A' * A;
disp(ATA);

% use gs to find the Gra
[Q, ~] = gs(A);
%display U
U=Q

%open a new figure called figure 3
figure(3);
%subplot each of 8 signals
%draw the 
for k=1:8
subplot(8,1,k);
stem(U(:,k)); axis([0 9 -0.5 0.5]); axis off; hold on;
end
for k=1:8
subplot(8,1,k);
plot(U(:,k));
end

%Compute the expansion coefficients ax with respect to the basis vectors of U.
a = U' * x;

%display the computed coefficients
disp('Expansion coefficients a:');
disp(a);

% Assuming 'a' is already defined

% Step 1: Find the absolute values of the coefficients in 'a'
abs_a = abs(a);

% Step 2: Identify the indices of the two largest entries in 'abs_a'
[~, indices] = maxk(abs_a, 2); % 'maxk' finds the k largest entries

% Step 3: Create a new vector 'a2' of the same size as 'a', with all zeros
a2 = zeros(size(a));

% Step 4: Assign the two largest entries from 'a' to 'a2'
a2(indices) = a(indices);

%Reconstruct x2 from a2
x2 = U * a2;

%plot x2 over Figure 1
figure(1); stem(x2,'b*'); plot(x2,'b');


% Identify the indices of the two largest entries in 'abs_a'
[~, indices] = maxk(abs_a, 4); % 'maxk' finds the k largest entries

% Create a new vector 'a4' of the same size as 'a', with all zeros
a4 = zeros(size(a));

% Assign the two largest entries from 'a' to 'a4'
a4(indices) = a(indices);

%Reconstruct x4 from a4
x4 = U * a4;

%plot x4 over Figure 1
figure(1); stem(x4,'k*'); plot(x4,'k');

%save the file as png
saveas(gcf, 'Figure_1.png')

%Reconstruct x8 from a
x8 = U * a;
% Compute relative errors
%error for x8
error_x8 = sqrt(sum((x - x8).^2) / sum(x.^2));
%error for x4
error_x4 = sqrt(sum((x - x4).^2) / sum(x.^2));
%error for x2
error_x2 = sqrt(sum((x - x2).^2) / sum(x.^2));

% printout the errors
fprintf('Relative error for x8: %f\n', error_x8);
fprintf('Relative error for x4: %f\n', error_x4);
fprintf('Relative error for x2: %f\n', error_x2);
%% Auxiliary function(s)
function [Q,R] = gs(A)
% Use CLASSICAL Gram Schmidt procedure to orthonormalizing the columns of 
% matrix A
% or equivalantly, perform the QR factorization of matrix A
% Input: matrix A with linearly independent columns
% Output: matrix Q, with orthonormal column vectors that spans span{col(A)}


    % n = # of columns of A 
    n = size(A,2);

    % Check if A has linearly independent columns. If not, end program and
    % return message 'Give a linearly independent set of vectors'
    % This part is optional

    if abs(det(transpose(A) * A)) < 10e-10
        fprintf('Give a linearly independent set of vectors \n') 
        return;
    end

    % Initialiize Q
    Q = A;
    R = zeros(n);

    % Compute q_1 by normalization
    R(1,1) = norm(Q(:,1));
    Q(:,1) = Q(:,1) / R(1,1);

    % Compute q_i by subtracting proj_{q_j}(A_i) from A_i, then normalize
    for i = 2:n
        
        % Substraction
        R(1:i-1,i) = Q(:,1:i-1)' * A(:,i);
        Q(:,i) = A(:,i) - Q(:,1:i-1) * R(1:i-1,i);
        
        % Normalization
        R(i,i) = norm(Q(:,i));
        Q(:,i) = Q(:,i) / R(i,i);
        
    end
end


