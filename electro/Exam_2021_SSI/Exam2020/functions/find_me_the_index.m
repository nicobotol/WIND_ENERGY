function index = find_me_the_index(vector, element)

    AA = vector == element;
    [~, index] = max(AA);

end