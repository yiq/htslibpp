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
                auto proxy {htsProxy(p)};
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
                auto proxy {HTSProxyIDPairContig(p)};
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
                auto proxy {htsProxy( p )};
                return proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::FORMAT) && (strcmp(proxy.key(), "AF") == 0);
            });

    // do I have the AO field?
    std::cout<<std::any_of(
            htsHeader<bcfHeader>::dictBegin(header, htsHeader<bcfHeader>::DictType::ID),
            htsHeader<bcfHeader>::dictEnd(header, htsHeader<bcfHeader>::DictType::ID),
            [](const auto& p) {
                auto proxy {htsProxy( p )};
                return proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::FORMAT) && (strcmp(proxy.key(), "AO") == 0);
            });

    return 0;
}
