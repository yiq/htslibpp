#include "htslibpp.h"
#include "htslibpp_proxies.h"
#include <iostream>
#include <algorithm>

using namespace YiCppLib::HTSLibpp;

int main(int argc, char* argv[]) {

    if(argc < 2) {
        std::cout<<"nice try. use a file name"<<std::endl;
        return 0;
    }

    auto fileHandle = htsOpen(argv[1], "r");
    if(not fileHandle) {
        std::cerr<<"unable to open file"<<std::endl;
        return 1;
    }

    auto header = htsHeader<bcfHeader>::read(fileHandle);

    auto hasDP = std::any_of(
            htsHeader<bcfHeader>::dictBegin(header, htsHeader<bcfHeader>::DictType::ID),
            htsHeader<bcfHeader>::dictEnd(header, htsHeader<bcfHeader>::DictType::ID),
            [](const auto& p) {
                auto proxy {htsProxy( p )};
                return proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::INFO) && (strcmp(proxy.key(), "DP") == 0);
            });
    std::cout<<hasDP<<std::endl;

    if(not hasDP) {
        std::cerr<<"The vcf file has no DP field in its INFO dictionary"<<std::endl;
        return 2;
    }

    int32_t *dp = (int32_t *)calloc(1, sizeof(int32_t));
    std::for_each(
            htsReader<bcfRecord>::begin(fileHandle, header),
            htsReader<bcfRecord>::end(fileHandle, header), 
            [&header, &dp](const auto& p) {
                bcf_unpack(p.get(), BCF_UN_INFO);
                int dp_size = 0;
                bcf_get_info_int32(header.get(), p.get(), "DP", &dp, &dp_size);
                std::cout<<dp[0]<<std::endl;
            });
    free(dp);

    return 0;
}





