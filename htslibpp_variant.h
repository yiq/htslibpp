// YiCppLib::HTSLibpp::BCF
//
// This file contains the c++14 wrapper of the VCF/BCF functionalities
// in htslib

#include "htslibpp.h"

#ifndef YICPPLIB_HTSLIBPP_BCF
#define YICPPLIB_HTSLIBPP_BCF

namespace YiCppLib {
    namespace HTSLibpp {
        // The next thing we need to do is to look at the content of the 
        // file once it is open. hts files are separated into two sections,
        // the header and the actual records. in htslib, the header is
        // represented by the type bcf_hdr_t. You get a pointer to this
        // struct by calling bcf_hdr_read, but you need to manually free
        // the struct with bcf_hdr_destroy. Using a smart pointer wrapper,
        // we can use RAII to manage the life cycle
        //using bcfHeader = HTS_UPTR(::bcf_hdr_t, void, bcf_hdr_destroy);
        using bcfHeader = HTS_UPTR(::bcf_hdr_t, bcf_hdr_destroy);

        // bcfHeader template generalization of HTSLibpp::htsHeader<T>
        // This struct groups functions that acts upon bcfHeader
        template<> struct htsHeader<bcfHeader> {

            inline static auto read(const htsFile& file) noexcept {
                return bcfHeader{ bcf_hdr_read(file.get()) };
            }

            // with a header, one can iterate over all its lines, which can
            // be of several types, currently there are 6 different types
            //   * FILTER,  e.g. ##FILTER=<ID=LowQual,Description="Low quality">
            //   * INFO,    e.g. ##INFO=<ID=AC_FIN,Number=A,Type=Integer,Description="Finnish Allele Counts">
            //   * FORMAT,  e.g. ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            //   * CONTIG,  e.g. ##contig=<ID=1,length=249250621>
            //   * STRUCT,  e.g. ##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
            //   * GENERAL, e.g. ##fileformat=VCFv4.1
            enum class LineType { FILTER, INFO, FORMAT, CONTIG, STRUCT, GENERAL };

            // header lines also contain type information, which can be one of
            //   * FLAG, essentially a boolean value with which presence or absence dictates some states to be true or false
            //   * INT,  e.g. DP=54
            //   * REAL, e.g. AF=0.48
            //   * STRING, e.g. culprit=MQ
            enum class DataType { FLAG, INT, REAL, STRING };

            // each of the header fields can have a different type of variable
            // length in each of the vcf records. 
            enum class VarLength { FIXED, VARIABLE, A, G, R };

            // Inside bcf_hdr_t, the records are stored in an array of bcf_hrec_t pointers, **hrec.
            // The size of the array hrec is specified by int nhrec. To iterate over all the records,
            // one need to perform pointer arithmatics with the type bcf_hrec_t**, but the actual
            // struct needs to be obtained by dereferencing twice of the type bcf_hrec_t**. 
            static inline auto begin(bcfHeader& hdr) { return hdr->hrec; }
            static inline auto end(bcfHeader& hdr) { return hdr->hrec + hdr->nhrec; }
            static inline const auto cbegin(bcfHeader& hdr) { return begin(hdr); }
            static inline const auto cend(bcfHeader& hdr) { return end(hdr); }

            // The headers are also organized into three different categories, which are
            //   * ID,      e.g. ##FILTER=<ID=LowQual,Description="Low quality">
            //   * CTG,     e.g. ##contig=<ID=1,length=249250621>
            //   * SAMPLE,  e.g. #CHROM ... FORMAT Sample1 Sample2
            enum class DictType { ID, CONTIG, SAMPLE };

            // One can choose to iterate over a particular dictionary
            static inline auto dictBegin(bcfHeader& hdr, DictType type) {
                switch(type) {
                    case DictType::ID:
                        return hdr->id[BCF_DT_ID];
                    case DictType::CONTIG:
                        return hdr->id[BCF_DT_CTG];
                    case DictType::SAMPLE:
                        return hdr->id[BCF_DT_SAMPLE];
                    default:
                        return static_cast<bcf_idpair_t *>(nullptr);
                }
            }

            static inline auto dictEnd(bcfHeader& hdr, DictType type) {
                switch(type) {
                    case DictType::ID:
                        return hdr->id[BCF_DT_ID] + hdr->n[BCF_DT_ID];
                    case DictType::CONTIG:
                        return hdr->id[BCF_DT_CTG] + hdr->n[BCF_DT_CTG];
                    case DictType::SAMPLE:
                        return hdr->id[BCF_DT_SAMPLE] + hdr->n[BCF_DT_SAMPLE];
                    default:
                        return static_cast<bcf_idpair_t *>(nullptr);
                }
            }
        };
   
        // Now let's get to actually reading records from a bcf / vcf file.
        // Each record is represented as bcf1_t struct in htslib. For the
        // purpose of performance, this struct is now allocated everytime
        // a record is fetched, but rather an already existing record,
        // being freshly allocated or not, is provided to the read function.
        using bcfRecord = HTS_UPTR(::bcf1_t, bcf_destroy);

        template<> struct htsReader<bcfRecord> {
            // Read the next bcf record form the file.
            static inline void read(htsFile& fp, const bcfHeader& hdr, bcfRecord& rec) {
                auto retVal = bcf_read(fp.get(), hdr.get(), rec.get());
                if(retVal == -1) rec.reset(nullptr);
            }

            // If one do not already possess a bcf1_t object, but simply want
            // the next bcf record, here is the approprate function
            static inline auto read(htsFile& fp, const bcfHeader& hdr) {
                bcfRecord rec{bcf_init()};
                read(fp, hdr, rec);
                return rec;
            }

            struct iterator : public std::iterator<std::input_iterator_tag, bcfRecord> {
                protected:
                    htsFile& fp;
                    const bcfHeader& hdr;
                    bcfRecord rec;

                public:
                    iterator(htsFile& fp, const bcfHeader& hdr) : fp(fp), hdr(hdr), rec(bcf_init()) {}
                    iterator(htsFile& fp, const bcfHeader& hdr, bcfRecord rec): fp(fp), hdr(hdr), rec(std::move(rec)) {}
                    virtual bcfRecord& operator*() { return rec; }

                    iterator& operator++()               { read(fp, hdr, rec); return *this; }
                    void operator++(int)                 { ++(*this); }

                    bool operator==(const iterator& rhs) { return rec.get() == nullptr && rhs.rec.get() == nullptr ; }
                    bool operator!=(const iterator& rhs) { return !(*this == rhs); }
            };

            static iterator begin(htsFile& fp, const bcfHeader& hdr) { return iterator{fp, hdr}; }
            static iterator end(htsFile& fp, const bcfHeader& hdr)   { return iterator{fp, hdr, std::move(bcfRecord{nullptr})}; }
        };
    }
}

namespace std {
    // iterator helper functions for bcfHeader
    auto inline begin(YiCppLib::HTSLibpp::bcfHeader& hdr) { return YiCppLib::HTSLibpp::htsHeader<YiCppLib::HTSLibpp::bcfHeader>::begin(hdr); }
    auto inline end(YiCppLib::HTSLibpp::bcfHeader& hdr) { return YiCppLib::HTSLibpp::htsHeader<YiCppLib::HTSLibpp::bcfHeader>::end(hdr); }
    auto inline cbegin(YiCppLib::HTSLibpp::bcfHeader& hdr) { return YiCppLib::HTSLibpp::htsHeader<YiCppLib::HTSLibpp::bcfHeader>::cbegin(hdr); }
    auto inline cend(YiCppLib::HTSLibpp::bcfHeader& hdr) { return YiCppLib::HTSLibpp::htsHeader<YiCppLib::HTSLibpp::bcfHeader>::cend(hdr); }

    // iterator helper functions for htsFile
    auto inline begin(YiCppLib::HTSLibpp::htsFile& fp, const YiCppLib::HTSLibpp::bcfHeader& hdr) { return YiCppLib::HTSLibpp::htsReader<YiCppLib::HTSLibpp::bcfRecord>::begin(fp, hdr);}
    auto inline end(YiCppLib::HTSLibpp::htsFile& fp, const YiCppLib::HTSLibpp::bcfHeader& hdr) { return YiCppLib::HTSLibpp::htsReader<YiCppLib::HTSLibpp::bcfRecord>::end(fp, hdr);}
}


// proxy classes
namespace YiCppLib {
    namespace HTSLibpp{

        /* The proxy around bcf_idpair_t, which normally contains fields defined
         * in FILTER / INFO / FORMAT header lines, but can also contain infomation
         * regarding CONTIG lines (see HTSProxyIDPairContig)
         */
        template<> struct HTSProxy<const bcf_idpair_t &> {
            const bcf_idpair_t& m_actual;

            HTSProxy(const bcf_idpair_t& actual) : m_actual(actual) {}
            inline auto key() const { return m_actual.key; }

            inline auto hasValueForLineType(htsHeader<bcfHeader>::LineType type) const {
                switch(type) {
                    case htsHeader<bcfHeader>::LineType::FILTER:
                        return m_actual.val->hrec[BCF_HL_FLT] != nullptr;
                    case htsHeader<bcfHeader>::LineType::INFO:
                        return m_actual.val->hrec[BCF_HL_INFO] != nullptr;
                    case htsHeader<bcfHeader>::LineType::FORMAT:
                        return m_actual.val->hrec[BCF_HL_FMT] != nullptr;
                    default:
                        return false;
                }
            }
        };

        /* If it is known that the bcf_idpair_t is referring to an contig
         * header line, this proxy should be used in instead. However such
         * choice cannot be based on the type of the raw object, but rather
         * made by the caller
         */
        struct HTSProxyIDPairContig : HTSProxy<const bcf_idpair_t &> {

            HTSProxyIDPairContig(const bcf_idpair_t& actual) : HTSProxy<const bcf_idpair_t&>(actual) {}

            inline auto hasValueForLineType(htsHeader<bcfHeader>::LineType type) {
                return type == htsHeader<bcfHeader>::LineType::CONTIG;
            }

            inline auto contigSize() { return m_actual.val->info[0]; }
        };
    }
}
#endif
