
reach_ind = 1:2*K;
removed = true;
while removed
    removed = false;
    for s = reach_ind
        if all(A(s,reach_ind) <= 0) && all(B(s,reach_ind) <= 0)
            removed = true;
            break;
        end
    end
    if removed
        reach_ind(reach_ind==s) = [];
    end
end     

Ap = A(reach_ind,reach_ind);
Bp = B(reach_ind,reach_ind);
Kp = length(reach_ind);

bad_input_ind = [];
for i=1:Kp
    if ones(1,Kp)*(Ap(:,i)+Bp(:,i)) ~= 1
        bad_input_ind(end+1) = i;
    end
end

req_input_ind = [];
for i=1:Kp
    if all(Ap(:,i) <= 0)
        req_input_ind(end+1) = i;
    end
end

Kf = zeros(Kp,Kp);
for j = req_input_ind
    Kf(j,j) = 1;
end
Af = Ap + Bp * Kf;

input_ind = 1:Kp;
input_ind(req_input_ind) = -1;
input_ind(bad_input_ind) = -1;
input_ind(input_ind == -1) = [];
Bf = Bp(:,input_ind);
orig_input_ind = reach_ind(input_ind);

orig_off_ind = reach_ind(reach_ind <= K);
orig_on_ind = reach_ind(reach_ind >= K+1);
off_ind = [];
on_ind = [];
for j=1:Kp
    if any(reach_ind(j) == orig_off_ind)
        off_ind(end+1) = j;
    end
    if any(reach_ind(j) == orig_on_ind)
        on_ind(end+1) = j;
    end
end


m = length(input_ind);