#include "parser_unittest.cpp"
#include "kmers_unittest.cpp"
#include "greedy_unittest.cpp"
#include "greedy_unittest_ac.cpp"
#include "generalized_simplitigs_unittest.cpp"
#include "generalized_simplitigs_ac_unittest.cpp"

#include "gtest/gtest.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
