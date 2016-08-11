#include "htslibpp.h"
#include "htslibpp_proxies.h"
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <ctime>

using namespace YiCppLib::HTSLibpp;

int main(int argc, char* argv[]) {

    if(argc<2) {
        std::cout<<"usage "<<argv[0]<<" <sample-rate> [input_filename]"<<std::endl;
        std::cout<<std::endl;
        std::cout<<"sample-rate is a probability a record will be chosen"<<std::endl;
        std::cout<<"if input_filename is not given, data is read from stdin"<<std::endl;
        return 0;
    }

    auto fileHandle = htsOpen(argc < 3 ? "-" : argv[2], "r");
    if(not fileHandle) {
        std::cerr<<"unable to open file"<<std::endl;
        return 1;
    }

    auto outFileHandle = htsOpen("-", "w");

    auto header = htsHeader<bcfHeader>::read(fileHandle);
    bcf_hdr_write(outFileHandle.get(), header.get());

    srand(time(NULL));

    const auto chance = std::stof(argv[1]);

    std::for_each(
            htsReader<bcfRecord>::begin(fileHandle, header),
            htsReader<bcfRecord>::end(fileHandle, header), 
            [&outFileHandle, &header, &chance](const auto& p) {
                auto diceRoll = rand();
                double roll = static_cast<double>(diceRoll) / static_cast<double>(RAND_MAX);
                if(roll <= chance) bcf_write(outFileHandle.get(), header.get(), p.get());
            });
    
    return 0;
}

