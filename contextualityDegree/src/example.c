#include "quadrics.h"

int main()
{
    main_header();

    quantum_assignment qa = subspaces(3, 1);
    bool sol[BV_LIMIT_CUSTOM(3)] = {0};
    int deg = geometry_contextuality_degree_custom(&qa, false, false, false, INVALID_LINES_HEURISTIC_SOLVER, sol);
    print("estimated degree of W(5,2):%d\n", deg);

    if (deg == 63) {
        quantum_assignment hex = quantum_assignment_from_invalid_contexts(qa, sol, false);

        print("split Cayley hexagon:\n");
        print_quantum_assignment(&hex);
    }

}