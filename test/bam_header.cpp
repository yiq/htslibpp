#include <gmock/gmock.h>
#include "../htslibpp.h"
#include "../htslibpp_alignment.h"

#include <algorithm>

using namespace YiCppLib::HTSLibpp;

class BamHeader: public testing::Test {
    public:
        YiCppLib::HTSLibpp::htsFile htsFileHandler = htsOpen("datasets/brca2.na12878.bam", "r");
};


TEST_F(BamHeader, CanCreateBamHeaderObject) {
    auto header = htsHeader<bamHeader>::read(htsFileHandler);
    ASSERT_NE(header.get(), nullptr);
}

TEST_F(BamHeader, CanIterateOverRecordsSequencially) {
    auto header = htsHeader<bamHeader>::read(htsFileHandler);
    size_t line_ct = 0;
    size_t hd_ct = 0;
    size_t sq_ct = 0;
    size_t pg_ct = 0;
    size_t rg_ct = 0;
    size_t co_ct = 0;
    std::for_each(
        htsHeader<bamHeader>::cbegin_l(header),
        htsHeader<bamHeader>::cend_l(header),
        [&](const auto& p) {
            line_ct++;

            auto tag = p.substr(0, 3);

            if(tag == "@HD")        hd_ct++;
            else if(tag == "@SQ")   sq_ct++;
            else if(tag == "@PG")   pg_ct++;
            else if(tag == "@RG")   rg_ct++;
            else if(tag == "@CO")   co_ct++;
        }
    );

    ASSERT_EQ(hd_ct, 1);
    ASSERT_EQ(sq_ct, 86);
    ASSERT_EQ(pg_ct, 28);
    ASSERT_EQ(rg_ct, 1);
    ASSERT_EQ(line_ct, 116);
}
