#include "htslibpp.h"
#include "htslibpp_variant.h"
#include <iostream>
#include <algorithm>

using namespace YiCppLib::HTSLibpp;

int main(int argc, char* argv[]) {
    
    auto fileHandle = htsOpen(argc < 2 ? "-" : argv[1], "r");

    if(fileHandle) {
        std::cout<<"hts file open of the format:";
        const auto hts_format = &fileHandle->format;
        auto c_string_dtor = [](char *str) {free(str);};
        auto format_str = std::unique_ptr<char, decltype(c_string_dtor)>(hts_format_description(hts_format), c_string_dtor);
        std::cout<<format_str.get()<<std::endl;
    }

    auto header = htsHeader<bcfHeader>::read(fileHandle);
    if(header) {
        std::cout<<"got bcf header of version: "<<bcf_hdr_get_version(header.get())<<std::endl;
    }

    std::for_each(std::cbegin(header), std::cend(header), [](const auto& p) {
            std::cout<<"got head of type "<<p->type<<std::endl;
        });

    auto dict = htsHeader<bcfHeader>::DictType::ID;
    std::for_each(
            htsHeader<bcfHeader>::dictBegin(header, dict),
            htsHeader<bcfHeader>::dictEnd(header, dict),
            [dict=dict](const auto& p) {
                auto proxy = htsProxy(p);
                std::cout<<proxy.key()<<": ";
                std::cout<<"isFilter<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::FILTER)<<">, ";
                std::cout<<"isInfo<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::INFO)<<">, ";
                std::cout<<"isFormat<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::FORMAT)<<">, ";
                std::cout<<"isContig<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::CONTIG)<<">, ";
                std::cout<<"isStruct<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::STRUCT)<<">, ";
                std::cout<<"isGeneral<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::GENERAL)<<">, ";
                
                std::cout<<std::endl;
            });

    dict = htsHeader<bcfHeader>::DictType::CONTIG;
    std::for_each(
            htsHeader<bcfHeader>::dictBegin(header, dict),
            htsHeader<bcfHeader>::dictEnd(header, dict),
            [](const auto& p) {
                auto proxy = HTSProxyIDPairContig(p);
                std::cout<<proxy.key()<<":";
                std::cout<<proxy.contigSize()<<". ";
                std::cout<<"isFilter<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::FILTER)<<">, ";
                std::cout<<"isInfo<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::INFO)<<">, ";
                std::cout<<"isFormat<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::FORMAT)<<">, ";
                std::cout<<"isContig<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::CONTIG)<<">, ";
                std::cout<<"isStruct<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::STRUCT)<<">, ";
                std::cout<<"isGeneral<"<<proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::GENERAL)<<">, ";
                std::cout<<std::endl;
            });

    // do I have AF fields in each sample, or in another word, part of 'FORMAT' line
    std::cout<<std::any_of(
            htsHeader<bcfHeader>::dictBegin(header, htsHeader<bcfHeader>::DictType::ID),
            htsHeader<bcfHeader>::dictEnd(header, htsHeader<bcfHeader>::DictType::ID),
            [](const auto& p) {
                auto proxy = htsProxy( p );
                return proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::FORMAT) && (strcmp(proxy.key(), "AF") == 0);
            });

    // do I have the AO field?
    std::cout<<std::any_of(
            htsHeader<bcfHeader>::dictBegin(header, htsHeader<bcfHeader>::DictType::ID),
            htsHeader<bcfHeader>::dictEnd(header, htsHeader<bcfHeader>::DictType::ID),
            [](const auto& p) {
                auto proxy = htsProxy( p );
                return proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::FORMAT) && (strcmp(proxy.key(), "AO") == 0);
            });

    //auto record = htsReader<bcfRecord>::read(fileHandle, header);
    //while(record.get() != nullptr) {
    //    std::cout<<record->rid<<"\t"<<record->pos<<std::endl;
    //    htsReader<bcfRecord>::read(fileHandle, header, record);
    //}

    int32_t *dp = (int32_t *)calloc(1, sizeof(int32_t));
    std::for_each(
            htsReader<bcfRecord>::begin(fileHandle, header),
            htsReader<bcfRecord>::end(fileHandle, header), 
            [&header, &dp](const auto& p) {
                bcf_unpack(p.get(), BCF_UN_INFO);
                int dp_size = 0;
                bcf_get_info_int32(header.get(), p.get(), "DP", &dp, &dp_size);
                std::cout<<p->rid<<"\t"<<p->pos<<"\t"<<p->unpacked<<"\tDP="<<dp[0]<<"DP_SIZE="<<dp_size<<std::endl;
            });
    free(dp);
        
    return 0;
}
