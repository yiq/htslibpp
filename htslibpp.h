#include <memory>
#include <string>
#include <functional>
#include <htslib/hts.h>
#include <htslib/vcf.h>

namespace YiCppLib {
    // a simple std::unique_ptr wrapper around traditional c-style pointers that
    // requires calling a function to release the internal data structure. 
    template<class T, class D> auto make_cptr_wrapper(T* p, D* d) { return std::unique_ptr<T, D>(p, d); }

    // a more complicated std::unique_ptr wrapper that allows destructor to be
    // defined at the time type alias is defined. Due to the syntax of C++
    // template system, the type of the destructor and the desctructor needs
    // to be specified individually, making a macro necessary to reduce typing
    template<class T, class D, D * d> struct _uptr_deleter { auto operator()(T* p) { if(p) d(p); } };
    template<class T, class D, D * d> using _hts_uptr = std::unique_ptr<T, _uptr_deleter<T, D, d>>;

    // A construct that will come in handy is an iterator wrapper around traditional
    // c-style pointer + size dynamic array.
    template<class T, class SIZE_T> struct _ptr_array_iterator : std::iterator<std::forward_iterator_tag, T> {
        protected:
            T * m_head;
            const SIZE_T m_size;
            SIZE_T m_curPos;

        public:
            _ptr_array_iterator(T * head, const SIZE_T& size, const SIZE_T& curPos = 0) : m_head(head), m_size(size), m_curPos(curPos) {}
            virtual T& operator*() { return *(m_head + m_curPos); }

            _ptr_array_iterator& operator++()    { ++m_curPos; return *this; }
            _ptr_array_iterator  operator++(int) { auto retVal = *this; ++(*this); return retVal; }
           
            bool operator==(const _ptr_array_iterator& rhs) { return m_curPos == rhs.m_curPos; }
            bool operator!=(const _ptr_array_iterator& rhs) { return !(*this == rhs); }

            static _ptr_array_iterator begin(T * head, SIZE_T& size) { return _ptr_array_iterator{head, size}; }
            static _ptr_array_iterator end(T * head, SIZE_T& size)   { return _ptr_array_iterator{head, size, size}; }
    };

    namespace HTSLibpp {
        // undefined generic htsHeader struct placeholder. Specific header types
        // will implement template specifications.
        template<class T> struct htsHeader;
    }
}

// Now the good stuff
#define HTS_UPTR(type, dtor) _hts_uptr<type, decltype(dtor), dtor>;
namespace YiCppLib {
    namespace HTSLibpp {

        // First order of business when dealing with hts files is to open it.
        // To follow the principle of Resource Acquisition Is Initialization,
        // or RAII for short, we define an type alias HTSLibpp::htsFile to
        // be an unique pointer to the underlying struct htsFile *, which
        // will be closed automatically upon going out-of-scope by calling
        // hts_close
        using htsFile = HTS_UPTR(::htsFile, hts_close);

        inline auto htsOpen(const std::string& filename, const std::string& mode) {
            return htsFile{ hts_open(filename.c_str(), mode.c_str()) };
        }

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
            enum class HeaderLine { FILTER, INFO, FORMAT, CONTIG, STRUCT, GENERAL };

            // header lines also contain type information, which can be one of
            //   * FLAG, essentially a boolean value with which presence or absence dictates some states to be true or false
            //   * INT,  e.g. DP=54
            //   * REAL, e.g. AF=0.48
            //   * STRING, e.g. culprit=MQ
            enum class HeaderType { FLAG, INT, REAL, STRING };

            enum class VariableLength { FIXED, VARIABLE, A, G, R };

            // Inside bcf_hdr_t, the records are stored in an array of bcf_hrec_t pointers, **hrec.
            // The size of the array hrec is specified by int nhrec. To iterate over all the records,
            // one need to perform pointer arithmatics with the type bcf_hrec_t**, but the actual
            // struct needs to be obtained by dereferencing twice of the type bcf_hrec_t**. 
            static inline auto begin(bcfHeader& hdr) { return hdr->hrec; }
            static inline auto end(bcfHeader& hdr) { return hdr->hrec + hdr->nhrec; }
        };
    }
    auto inline begin(HTSLibpp::bcfHeader& hdr) { return HTSLibpp::htsHeader<HTSLibpp::bcfHeader>::begin(hdr); }
    auto inline end(HTSLibpp::bcfHeader& hdr) { return HTSLibpp::htsHeader<HTSLibpp::bcfHeader>::end(hdr); }

}
