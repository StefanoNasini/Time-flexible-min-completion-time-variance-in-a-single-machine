function y = fact(x)

if x == 0
    y = 1;
    return
end

y = x*fact(x-1);

end