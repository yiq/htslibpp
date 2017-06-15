#include <memory>
#include <string>
#include <functional>
#include <htslib/hts.h>
#include <htslib/vcf.h>

#ifndef YICPPLIB_HTSLIBPP_HTSLIBPP
#define YICPPLIB_HTSLIBPP_HTSLIBPP

template<class T> class TD;     // A crude instrument to check compiler deduced type

namespace YiCppLib {
    // a simple std::unique_ptr wrapper around traditional c-style pointers that
    // requires calling a function to release the internal data structure. 
    template<class PtrT, class DtorT> auto make_cptr_wrapper(PtrT* p, DtorT* d) { return std::unique_ptr<PtrT, DtorT>(p, d); }

    // a more complicated std::unique_ptr wrapper that allows destructor to be
    // defined at the time type alias is defined. Due to the syntax of C++
    // template system, the type of the destructor and the desctructor needs
    // to be specified individually, making a macro necessary to reduce typing
    template<class PtrT, class DtorT, DtorT * d> struct _uptr_deleter { auto operator()(PtrT* p) { if(p) d(p); } };
    template<class PtrT, class DtorT, DtorT * d> using _uptr_with_dtor = std::unique_ptr<PtrT, _uptr_deleter<PtrT, DtorT, d>>;

    // A construct that will come in handy is an iterator wrapper around traditional
    // c-style pointer + size dynamic array.
    template<class T, class SIZE_T> struct _ptr_array_iterator : std::iterator<std::bidirectional_iterator_tag, T> {
        protected:
            T * m_head;
            const SIZE_T m_size;
            SIZE_T m_curPos;

        public:
            _ptr_array_iterator(T * head, const SIZE_T& size, const SIZE_T& curPos = 0) : m_head(head), m_size(size), m_curPos(curPos) {}
            virtual T& operator*() { return *(m_head + m_curPos); }

            _ptr_array_iterator& operator++()    { ++m_curPos; return *this; }
            _ptr_array_iterator  operator++(int) { auto retVal = *this; ++(*this); return retVal; }

            _ptr_array_iterator& operator--()    { if(m_curPos > 0) --m_curPos; return *this; }
            _ptr_array_iterator  operator--(int) { auto retVal = *this; --(*this); return retVal; }
           
            bool operator==(const _ptr_array_iterator& rhs) { return m_curPos == rhs.m_curPos; }
            bool operator!=(const _ptr_array_iterator& rhs) { return !(*this == rhs); }

            static _ptr_array_iterator begin(T * head, SIZE_T& size) { return _ptr_array_iterator{head, size}; }
            static _ptr_array_iterator end(T * head, SIZE_T& size)   { return _ptr_array_iterator{head, size, size}; }
    };
    namespace HTSLibpp {
        // undefined generic htsHeader struct placeholder. Specific header types
        // will implement template specifications.
        template<class T> struct htsHeader;

        // undefined generic htsReader struct placeholder. Specific reader types
        // will implement template specifications.
        template<class T> struct htsReader;

        
        
        // proxy base templates

        /* The uninstanciated holder struct
         */
        template<class T> struct HTSProxy;

        /* The convenient function that makes proxy objects based on the type
         * of the raw object passed in as argument
         */
        template<class T>
        inline auto htsProxy(const T& actual) { return HTSProxy<decltype(actual)>(actual); }
    }

}

// Now the good stuff
#define HTS_UPTR(type, dtor) _uptr_with_dtor<type, decltype(dtor), dtor>;
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
    }
}

#endif
