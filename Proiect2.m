%% Identification part
clear all;
load('iddata-10.mat');
%extracting the identification values values
ts = id.Ts;
u_id = id.u;
y_id = id.y;
% checking received values
% subplot(211);
% plot(u_id);title("uid");
% subplot(212);
% plot(y_id); title("yid");

%getting values for na=nb from keyboard
na = input('Enter na and nb: ');
nb = na;
nk = 1;
%obtaining the number of elements
M=length(y_id);
phi=ones(M,na+nb);
%phi matrix creation
% for y(k-na)
for i=1:M
    for j=1:na
        if((i-j)<=0)
            phi(i,j)=0;
        elseif((i-j)>0) 
            phi(i,j)=-y_id(i-j);
        end
    end
end
% for y(k-nb)
for i=1:M
    for j=na+1:nb+na
        if((i+na-j)<=0)
            phi(i,j)=0;
        elseif((i+na-j)>0) 
            phi(i,j)=u_id(i+na-j);
        end
    end
end
%initial theta and yhat
theta = phi\y_id;
yhat = phi * theta;
%getting m value from keyboard input
m = input('Enter m: ');
%inputing the phi values in a cell array
x = cell(1, na+nb);
for i = 1:na+nb 
    x{i} = phi(:,i);
end

phiIden = []; 
% creating a matrix with the polynomials up to degree m
for i = 1:M %going up to the length of M
    h = [];%initializating empty matrix h
    k = 1;%counter
    nIndex = 1;
    if m == 1%degree=1 special case
       phiIden = phi;
    else   
        while m >= k
              for k = 1:m
                  for l = 1:na+nb
                      h = [h, x{l}(i)^(k)]; %individual elements up to necessary degree rank 
                  end
              end    
              if k == m    
                 for f = 1:na+nb-1
                     for g = f+1:na+nb
                         h = [h, x{f}(i)*x{g}(i)]; % multiplication for 2 variables 
                     end
                 end
                 for k = 2:m-1 
                     for f = 1:na+nb
                         for g = 1:na+nb
                             if f == g 
                                continue
                             end
                             for v = 1:k - 1 
                                 h = [h, (x{f}(i))^(k)*(x{g}(i))^(v)]; % multiplication for 2 variables and higher degree rank than 1
                             end
                         end
                     end  
                     for k = 2:m-1
                         for f = 1:na+nb-2
                             for g = f+1:na+nb-1
                                 for s = g+1:na+nb
                                    h = [h, x{f}(i)*x{g}(i)*x{s}(i)]; %multiplication for multiple variables
                                 end
                             end
                         end
                         if k == m - 1
                            k = m + 1;%modiyfying condition
                         end
                     end 
                 end
              end 
        end  
        phiIden(end+1, 1:length(h)) = h; %phi matrix for identification
    end 
end

theta2 = phiIden\yhat; %obtaining larger theta using left matrix division

yhatIden = phiIden*theta2; %final yhat

e = [];%error values array
for i = 1: length(yhatIden)
    e = [e, y_id(i) - yhatIden(i)];
end
mse = mean(e.^2); %mean square error value
figure(1)
plot(yhatIden);hold
plot(y_id);title("Identification Plot");
%% Validation
%for validation we use the same steps as in the identification phase but using the validation values
clear all;
load('iddata-10.mat');
ts = id.Ts;
u_val = val.u;
y_val = val.y;

% subplot(211);
% plot(u_id);title("uid");
% subplot(212);
% plot(y_id); title("yid");
na = input('Enter na and nb: ');
nb = na;
nk = 1;
M2=length(y_val);
phi=ones(M2,na+nb);

% for y(k-na)
for i=1:M2
    for j=1:na
        if((i-j)<=0)
            phi(i,j)=0;
        elseif((i-j)>0) 
            phi(i,j)=-y_val(i-j);
        end
    end
end
% for y(k-nb)
for i=1:M2
    for j=na+1:nb+na
        if((i+na-j)<=0)
            phi(i,j)=0;
        elseif((i+na-j)>0) 
            phi(i,j)=u_val(i+na-j);
        end
    end
end

theta = phi\y_val;
yhat = phi * theta;

m = input('Enter m: ');

x = cell(1, na+nb);
for i = 1:na+nb 
    x{i} = phi(:,i);
end

phiVal = []; 
for i = 1:M2 
    h = [];
    k = 1;
    nIndex = 1;
    if m == 1
       phiVal = phi;
    else   
        while m >= k
              for k = 1:m
                  for l = 1:na+nb
                      h = [h, x{l}(i)^(k)];
                  end
              end    
              if k == m    
                 for f = 1:na+nb-1
                     for g = f+1:na+nb
                         h = [h, x{f}(i)*x{g}(i)];
                     end
                 end
                 for k = 2:m-1 
                     for f = 1:na+nb
                         for g = 1:na+nb
                             if f == g 
                                continue
                             end
                             for v = 1:k - 1 
                                 h = [h, (x{f}(i))^(k)*(x{g}(i))^(v)];
                             end
                         end
                     end  
                     for k = 2:m-1
                         for f = 1:na+nb-2
                             for g = f+1:na+nb-1
                                 for s = g+1:na+nb
                                    h = [h, x{f}(i)*x{g}(i)*x{s}(i)];
                                 end
                             end
                         end
                         if k == m - 1
                            k = m + 1;
                         end
                     end 
                 end
              end 
        end  
        phiVal(end+1, 1:length(h)) = h; 
    end 
end

theta2 = phiVal\yhat;

yhatVal = phiVal*theta2;

e = [];
for i = 1: length(yhatVal)
    e = [e, y_val(i) - yhatVal(i)];
end
mse = mean(e.^2);
figure(2)
plot(yhatVal);hold
plot(y_val);title("Validation Plot");