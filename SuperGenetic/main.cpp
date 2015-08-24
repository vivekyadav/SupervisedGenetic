//
//  main.cpp
//  SuperGenetic
//
//  Created by D, Vivek on 17/06/15.
//  Copyright (c) 2015 Zero. All rights reserved.
//

#include <iostream>
#include "SuperGenetic.cpp"

int main(int argc, const char * argv[]) {
    SuperGenetic::TSPSolver solver = SuperGenetic::TSPSolver(1000, 100, 1000, 5, 1);
    solver.solve();
	//std::vector<int> result = solver.nearest_neighbour_solver(std::vector<int> {5, 6, 2, 3, 9});
	//for (int i : result) {
	//	std::cout << i << " ";
	//}
	getchar();
    return 0;
}
