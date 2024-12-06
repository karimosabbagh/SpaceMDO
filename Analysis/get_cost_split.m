function cost_split = get_cost_split(output, PB)
    [results, evaluations] =  evaluate_subsystems(output.x, PB);
    cost_split = evaluations.cost_split;
end