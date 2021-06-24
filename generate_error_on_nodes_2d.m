function error = generate_error_on_nodes_2d(exact_solution, PbU, U)

%Evaluation of the given analytical solution at all the nodes

y = size(PbU,2);

for i = 1:y
    vector_analytic(i) = feval(exact_solution,PbU(1,i),PbU(2,i));
end

error_vector = U - vector_analytic;

error = max(abs(error_vector));