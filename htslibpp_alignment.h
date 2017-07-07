// YiCppLib::HTSLibpp::Alignment
//
// This file contains the c++14 wrapper of the SAM/BAM/CRAM functionalities
// in htslib

#include "htslibpp.h"
#include <htslib/sam.h>
#include <string.h>
#ifndef YICPPLIB_HTSLIBPP_ALIGNMENT
#define YICPPLIB_HTSLIBPP_ALIGNMENT

// A BAM/SAM/CRAM file is separated into two parts, the header, and the body.
// The header consists of metadata describe things like
//   * file version
//   * references
//   * programs
//   * etc...
//
// The body consists of actual alignments, and is represented by bamRecord 
// objects. This part can be iterated sequencially, either for the entire
// file or a specific region.


// --- BAM HEADER --- //

namespace YiCppLib {
    namespace HTSLibpp {
        // The smart pointer wrapper to the SAM/BAM header
        using bamHeader = HTS_UPTR(::bam_hdr_t, bam_hdr_destroy);

        // bamHeader template generalization of HTSLibpp::htsHeader<T>
        template<> struct htsHeader<bamHeader> {

            inline static auto read(const htsFile& file) noexcept {
                return bamHeader{ bam_hdr_read(file->fp.bgzf) };
            }

            struct line_iterator : public std::iterator<std::forward_iterator_tag, std::string> {
                protected:
                    size_t pos;
                    const char * text;
                    const size_t l_text;

                public:
                    line_iterator(const bamHeader& header): line_iterator(header, 0) {}
                    line_iterator(const bamHeader& header, size_t pos): pos(pos), text(header->text), l_text(header->l_text) {}
                    virtual std::string operator*() { 
                        size_t start_pos = pos;
                        size_t end_pos = strchr(text+start_pos, '\n') - text;
                        return std::string(text+start_pos, text+end_pos); 
                    }

                    iterator& operator++() { pos = (strchr(text + pos, '\n') - text)+ 1; return *this; }
                    void operator++(int) { ++(*this); }

                    bool operator==(const line_iterator& rhs) { return rhs.pos == pos; }
                    bool operator!=(const line_iterator& rhs) { return !(*this == rhs); }
            };

            inline static auto cbegin_l(const bamHeader& header) noexcept { return line_iterator(header); }
            inline static auto cend_l(const bamHeader& header) noexcept { return line_iterator(header, header->l_text); }

        };
    }
}


// --- BAM RECORDS --- //
namespace YiCppLib {
    namespace HTSLibpp {
        // Each alignment in HTSLib is represented as bam1_t struct. For the
        // purpose of performance, this struct is not allocated everytime
        // a record is fetched, but rather an already existing record,
        // being freshly allocated or not, is provided to the read function
        using bamRecord = HTS_UPTR(::bam1_t, bam_destroy1);

        template<> struct htsReader<bamRecord> {
            // read the next bam record from the file.
            static inline void read(htsFile& fp, bamRecord& rec) {
                auto retVal = bam_read1(fp->fp.bgzf, rec.get());
                if(retVal == -1) rec.reset(nullptr);
            }

            // if one does not already process a bamRecord object, but simply wants
            // the next bam record, here is the approprate function
            static inline auto read(htsFile& fp) {
                bamRecord rec{bam_init1()};
                read(fp, rec);
                return rec;
            }

            // sequencial iterator
            struct iterator : public std::iterator<std::input_iterator_tag, bamRecord> {
                protected:
                    htsFile& fp;
                    bamRecord rec;

                public:
                    iterator(htsFile& fp): fp(fp), rec(bam_init1()) {}
                    iterator(htsFile& fp, bamRecord rec): fp(fp), rec(std::move(rec)) {}
                    virtual bamRecord& operator*() { return rec; }

                    iterator& operator++()               { read(fp, rec); return *this; }
                    void operator++(int)                 { ++(*this); }

                    bool operator==(const iterator& rhs) { return rec.get() == nullptr && rhs.rec.get() == nullptr; }
                    bool operator!=(const iterator& rhs) { return !(*this == rhs); }
            };

            static iterator begin(htsFile& fp) { return iterator{fp}; }
            static iterator end(htsFile& fp) { return iterator{fp, std::move(bamRecord{nullptr})}; }
        };
    }
}

// proxy classes
namespace YiCppLib {
    namespace HTSLibpp {
        // the proxy around bamRecord, which describes an alignment.
        template<> struct HTSProxy<const bamRecord &> {
            const bamRecord& m_actual;

            HTSProxy(const bamRecord& actual): m_actual(actual) {}

            inline auto chrID() const     { return m_actual->core.tid;  }
            inline auto pos() const       { return m_actual->core.pos;  }
            inline auto qual() const      { return m_actual->core.qual; }
            inline auto mateChrID() const { return m_actual->core.mtid; }
            inline auto matePos() const   { return m_actual->core.mpos; }

            // more expensive extractions
            //inline auto queryName() const { return std::string((char*)m_actual->data); } 
            //inline auto sequence() const { return std::string((char*)(m_actual->data + (m_actual->core.n_cigar)<<2 + m_actual->core.l_qname) ); }
            //inline auto baseQualities() const { return std::string((char*)(m_actual->data + (m_actual->core.n_cigar)<<2 + m_actual->core.l_qname + (m_actual->core.l_qseq + 1) >> 1)); } 
            // TO-BE-IMPLEMENTED
            // inline auto auxiliary() const;
            // inline auto auxiliary_length() const;
            // inline auto get_base(uint32_t i) const;
        };
    }
}

#endif
