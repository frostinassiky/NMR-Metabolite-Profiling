function F = sec_order_fit(x,h_sam,x_sam)
    num = length(x)/2;
    A = zeros(num);
    B = zeros(num);
    h_star = x(1:num);
    x_star = x(num+1:end);
    for k = 1:num
        for l = 1:num
            A(k,l) = lorezf(x_sam(k),x_star(l),1);
            B(k,l) = A(k,l)^2 * (x_sam(k)-x_star(l));
        end
    end
    F1 = A * h_star - h_sam;
    F2 = B * h_star;
    F = [F1;F2];
end