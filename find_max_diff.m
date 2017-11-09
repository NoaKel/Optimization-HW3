function [ max_diff_index ] = find_max_diff( vector )

n=length(vector);
diff=vector(1:n-1)-(vector(2:n));
[~,max_diff_index]=max(abs(diff));

end

