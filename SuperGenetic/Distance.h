#pragma once

#include<vector>
#include<math.h>

namespace Distance {
    template <typename T>
    inline double euclidian_distance(const std::vector<T>& point1, const std::vector<T>& point2) {
        check_dimension_equality(point1, point2);
        double distance = 0;
        for (int i = 0; i < point1.size(); ++i) {
            distance += pow(point1.at(i) - point2.at(i), 2);
        }
        return sqrt(distance);
    }

    template <typename T>
    inline int hamming_distance(const std::vector<T>& point1, const std::vector<T>& point2) {
        check_dimension_equality(point1, point2);
        int distance = 0;
        for (int i = 0; i < point1.size(); ++i) {
            if (point1.at(i) != point2.at(i)) {
                ++distance;
            }
        }
        return distance;
    }

    template <typename T>
    inline void check_dimension_equality(const std::vector<T>& point1, const std::vector<T>& point2) {
        if (point1.size() != point2.size()) {
            throw std::exception("Given points are of unequal dimension");
        }
    }
}
