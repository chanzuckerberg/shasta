#ifndef CZI_NANOPORE2_TIME_STAMP_HPP
#define CZI_NANOPORE2_TIME_STAMP_HPP

#include "boost/date_time/posix_time/posix_time.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {
        inline std::ostream& timestamp(std::ostream& s)
        {
            s << boost::posix_time::microsec_clock::local_time() << " ";
            return s;
        }
    }
}


#endif

