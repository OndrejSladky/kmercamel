#include "parser_unittest.h"
#include "kmers_unittest.h"
#include "global_unittest.h"
#include "global_ac_unittest.h"
#include "local_unittest.h"
#include "local_ac_unittest.h"
#include "streaming_unittest.h"

#include "gtest/gtest.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
