function abs_LD = abs_LD(x,M)
% % M should be row vector of size 1xk
 column_count = size(M,2);
 first_non_zero = 0; 
 if ne(sum(x), 0)
  first_non_zero = sign(x); 
 else    
  for k=1:column_count
       S = sum(M(:,k));
       if ne(S,0)
          first_non_zero =  sign(M(:,k));
          break
       end
  end
 end
 fsign = first_non_zero';
 abs_LD = fsign*M;
end