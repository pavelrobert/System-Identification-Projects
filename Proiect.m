load('C:\Users\Roby\Desktop\Uni Stuff\An3Sem1\SystemId\Proiect1\proj_fit_02.mat')
%mesh(id.X{1,1}, id.X{2,1}, id.Y);
%% initializing the identification values
x1_id = zeros(1,41);
x1_id = id.X{1,1};
x2_id = zeros(1,41);
x2_id = id.X{2,1};
yId = zeros(41,41);
yId = id.Y;

%% making the phi matrix as a cell array of vectors 
phi = cell(41, 41);
phi2 = []; % initializing the matrix that will be used for calculating the theta

m = input('Enter m: '); % user inputs the "m" variable

for i = 1:length(x1_id)     
    for j = 1:length(x1_id)     % 2 FORs for making the phi matrix
        h = [];                 % initializing a vector blank "h"
        k = 1;                  % initializing the variable "k", wich will be used for the powers of the Xs
        nIndex = 0;             % the "nIndex" variable is used for using the next index of the blank vector's last inputed element
        while m >= k            
            if nIndex == 0
               h = [h, 1];      % putting the number 1 at the start of each vector
            end
            nIndex = length(h) + 1;     
            if mod(nIndex, 2) == 0      % if the index of the vector is even we will use numbers from the x1_id vector, x2_id otherwise
               h = [h, x1_id(i)^k];
            else
               h = [h, x2_id(j)^k];
               k = k + 1;               % after each x2 value has been written in the h vector the polynomial 
                                        % rank increases as shown in the project description
            end
            if k == m + 1               % after finishing writing the powers of x1 and x2, the multiplications start
               h = [h, x1_id(i)*x2_id(j)];                      % the first multiplication done will always be x1*x2
               for k = 2:m-1                                    % we use a for loop to multiply the maximum k while the v and w are used for 
                   for v = 1:k - 1                              % ascending polynomial ranks and descending ones for example:
                       h = [h, x1_id(i)^(k)*x2_id(j)^(v)];      % x1^3*x2   x1^3*x2^2 x1*x2^3 x1^2*x2^3
                   end
                   for w = 1:k - 1   
                       h = [h, x1_id(i)^(w)*x2_id(j)^(k)];
                   end
                   if k == m - 1                                % if used for exiting the while loop 
                       k = m + 1;
                   end    
               end
            end 
        end
        phi{i,j} = h;                                           % adding the h vector into the phi cell array
        phi2(end+1, 1:length(h)) = h;                           
    end  
end

yId2 = reshape(yId, 1, []);
yId2 = yId2';
theta = phi2\yId2;

yhat = phi2*theta;
yhat_rs = reshape(yhat, 41, []);

mesh(id.X{1,1}, id.X{2,1}, yId);
figure
mesh(id.X{1,1}, id.X{2,1}, yhat_rs);

epsilon = zeros(41, 41);
for i = 1:41
    for j = 1:41
        epsilon(i,j) = abs((yId(i,j) - yhat_rs(i,j))^2);
    end
end

MSE_id = mean(epsilon, 'all');

%%
x1_val = zeros(1,71);
x1_val = val.X{1,1};
x2_val = zeros(1,71);
x2_val = val.X{2,1};
yVal = zeros(71,71);
yVal = val.Y;

phi_val = cell(71, 71);
phi2_val = [];

m = input('Enter m: ');

for i = 1:length(x1_val)
    for j = 1:length(x1_val)
        h = [];
        k = 1;
        nIndex = 0;
        while m >= k
            if nIndex == 0
               h = [h, 1];
            end
            nIndex = length(h) + 1;
            if mod(nIndex, 2) == 0
               h = [h, x1_val(i)^k];
            else
               h = [h, x2_val(j)^k];
               k = k + 1;
            end
            if k == m + 1
               h = [h, x1_val(i)*x2_val(j)];
               for k = 2:m-1
                   for v = 1:k - 1
                       h = [h, x1_val(i)^(k)*x2_val(j)^(v)];
                   end
                   for w = 1:k - 1   
                       h = [h, x1_val(i)^(w)*x2_val(j)^(k)];
                   end
                   if k == m - 1
                       k = m + 1;
                   end    
               end
            end 
        end
        phi_val{i,j} = h;
        phi2_val(end+1, 1:length(h)) = h;
    end  
end

yVal2 = reshape(yVal, 1, []);
yVal2 = yVal2';
theta2 = phi2_val\yVal2;

yhat2 = phi2_val*theta2;
yhat2_rs = reshape(yhat2, 71, []);


mesh(val.X{1,1}, val.X{2,1}, yVal);
hold
surf(val.X{1,1}, val.X{2,1}, yhat2_rs);

epsilon2 = zeros(71, 71);
for i = 1:71
    for j = 1:71
        epsilon2(i,j) = abs((yVal(i,j) - yhat2_rs(i,j))^2);
    end
end

MSE_val = mean(epsilon2, 'all');
