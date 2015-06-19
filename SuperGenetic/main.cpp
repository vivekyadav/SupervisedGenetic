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
    SuperGenetic::TSPSolver solver = SuperGenetic::TSPSolver(14, 10, 100, 5, 1);
    solver.solve();
    return 0;
}
