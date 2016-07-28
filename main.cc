#include "htslibpp.h"
#include <iostream>

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

    for(auto& hrec : header) std::cout<<hrec->type<<std::endl;

    return 0;

}
