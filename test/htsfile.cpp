#include<gmock/gmock.h>
#include "../htslibpp.h"

using namespace YiCppLib::HTSLibpp;

TEST(HTSLibpp, CanCreateHTSFileObjectWithBamFile) {
    auto htsFileHandle = htsOpen("datasets/brca2.na12878.bam", "r");
    ASSERT_NE(htsFileHandle.get(), nullptr);
}
