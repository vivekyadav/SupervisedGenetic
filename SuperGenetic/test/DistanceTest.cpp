#include "gtest/gtest.h"
#include "Distance.h"

#include <vector>

class DistanceTest : public testing::Test {
protected:
    virtual void setUp() {
    }

    const std::vector<double> point1 = { 1, 2 };
    const std::vector<double> point2 = { 4, 6 };
};

TEST_F(DistanceTest, euclidian_distance_test_int_type) {
    std::vector<int> point1 = { 1, 2 };
    std::vector<int> point2 = { 4, 6 };
    ASSERT_EQ(5, Distance::euclidian_distance<int>(point1, point2));
}

TEST_F(DistanceTest, euclidian_distance_test_double_type) {
    ASSERT_EQ(5, Distance::euclidian_distance<double>(point1, point2));
}

TEST_F(DistanceTest, euclidian_distance_unequal_dimension) {
    std::vector<double> point_with_higher_dimension = { 4, 6, 5 };
    ASSERT_THROW(Distance::euclidian_distance<double>(point1, point_with_higher_dimension), std::exception);
}

TEST_F(DistanceTest, hamming_distance_test_double_type) {
    ASSERT_EQ(2, Distance::hamming_distance<double>(point1, point2));
}

TEST_F(DistanceTest, hamming_distance_test_identical_points) {
    ASSERT_EQ(0, Distance::hamming_distance<double>(point1, point1));
}

TEST_F(DistanceTest, pseudo_euclidian_distance_double_type) {
    ASSERT_EQ(1, Distance::pseudo_euclidian_distance(point1, point2));
}

