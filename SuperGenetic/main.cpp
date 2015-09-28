//
//  main.cpp
//  SuperGenetic
//
//  Created by D, Vivek on 17/06/15.
//  Copyright (c) 2015 Zero. All rights reserved.
//

#include <iostream>
#include "SuperGenetic.cpp"
#include "tsplib.h"

int main(int argc, const char * argv[]) {
	TspLib::TspDataLoader data_loader;
	TspLib::TspData data = data_loader.load("C:\\Users\\vivek\\Desktop\\TSP\\att48.tsp");
    SuperGenetic::TSPSolver solver = SuperGenetic::TSPSolver(data.get_all_coordinates(), data.get_cities().size(), 5, 200, 10, 3);
    solver.solve();
	//std::vector<int> result = solver.nearest_neighbour_solver(std::vector<int> {5, 6, 2, 3, 9});
	//for (int i : result) {
	//	std::cout << i << " ";
	//}
	getchar();
    return 0;
}
