#include "parser_unittest.h"
#include "kmers_unittest.h"
#include "greedy_unittest.h"
#include "greedy_unittest_ac.h"
#include "generalized_simplitigs_unittest.h"
#include "generalized_simplitigs_ac_unittest.h"

#include "gtest/gtest.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
