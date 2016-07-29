/* @file htslibpp_proxies.h
 *
 * This file defines several proxy objects that wraps around raw hts types.
 * These proxy objects usually will hold a constant reference to the under-
 * lying data structure, and forward calls using hts standard functions
 */

#ifndef YICPPLIB_HTSLIBPP_PROXIES_H
#define YICPPLIB_HTSLIBPP_PROXIES_H

namespace YiCppLib {
    namespace HTSLibpp{

        /* The uninstanciated holder struct
         */
        template<class T> struct HTSProxy;

        /* The convenient function that makes proxy objects based on the type
         * of the raw object passed in as argument
         */
        template<class T>
        inline auto htsProxy(const T& actual) { return HTSProxy<decltype(actual)>(actual); }

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
