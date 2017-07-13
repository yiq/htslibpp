#include <gmock/gmock.h>
#include "../htslibpp.h"
#include "../htslibpp_alignment.h"

#include <algorithm>

using namespace YiCppLib::HTSLibpp;

class BamRecord : public testing::Test {
    public:
        YiCppLib::HTSLibpp::htsFile htsFileHandler = htsOpen("datasets/brca2.na12878.bam", "r");
};

TEST_F(BamRecord, CanIterateRecordsSequencially) {
    size_t read_count = 0;
    auto header = htsHeader<bamHeader>::read(htsFileHandler);

    std::for_each(
            htsReader<bamRecord>::begin(htsFileHandler, header),
            htsReader<bamRecord>::end(htsFileHandler, header),
            [&read_count](auto& r) { read_count++; });

    ASSERT_EQ(read_count, 45256);
}

TEST_F(BamRecord, CanCountRecords) {
    auto header = htsHeader<bamHeader>::read(htsFileHandler);

    auto read_count = std::distance(
            htsReader<bamRecord>::begin(htsFileHandler, header),
            htsReader<bamRecord>::end(htsFileHandler, header));

    ASSERT_EQ(read_count, 45256);
}

TEST_F(BamRecord, CanIterateRecordUsingRangeExpression) {
    size_t read_count = 0;
    auto header = htsHeader<bamHeader>::read(htsFileHandler);

    for(auto &r : htsReader<bamRecord>::range(htsFileHandler, header)) read_count++;

    ASSERT_EQ(read_count, 45256);
}

TEST_F(BamRecord, CanIterateRegionSequencially) {
    size_t read_count = 0;
    auto header = htsHeader<bamHeader>::read(htsFileHandler);
    auto index  = htsIndexOpen("datasets/brca2.na12878.bam", "datasets/brca2.na12878.bam.bai");
    const std::string region = "13:32900000-32950000";

    std::for_each(
            htsReader<bamRecord>::begin(htsFileHandler, header, index, region),
            htsReader<bamRecord>::end(htsFileHandler, header, index, region),
            [&read_count](const auto& p) { read_count++; });

    ASSERT_EQ(read_count, 27112);
}

//TEST_F(BamRecord, CanCountRecordsInRegion) {
//    auto header = htsHeader<bamHeader>::read(htsFileHandler);
//    auto index  = htsIndexOpen("datasets/brca2.na12878.bam", "datasets/brca2.na12878.bam.bai");
//    const std::string region = "13:32900000-32950000";
//    auto read_count = std::distance(
//            htsReader<bamRecord>::begin(htsFileHandler, header, index, region),
//            htsReader<bamRecord>::end(htsFileHandler, header, index, region));
//    ASSERT_EQ(read_count, 27112);
//}

TEST_F(BamRecord, CanIterateRegionUsingRangeExpression) {
    size_t read_count = 0;
    auto header = htsHeader<bamHeader>::read(htsFileHandler);
    auto index  = htsIndexOpen("datasets/brca2.na12878.bam", "datasets/brca2.na12878.bam.bai");
    const std::string region = "13:32900000-32950000";

    for(auto &r : htsReader<bamRecord>::range(htsFileHandler, header, index, region)) read_count++;

    ASSERT_EQ(read_count, 27112);
}

TEST_F(BamRecord, CanGetQueryName) {

    auto header = htsHeader<bamHeader>::read(htsFileHandler);
    auto firstRead = htsReader<bamRecord>::begin(htsFileHandler, header);
    auto proxy = htsProxy(*firstRead);

    ASSERT_EQ(proxy.queryName(), "ERR194147.537888192");
    //ASSERT_EQ(proxy.cigar(), "101M");
    
}
