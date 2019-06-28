#ifndef SHASTA_TIME_STAMP_HPP
#define SHASTA_TIME_STAMP_HPP

#include "boost/date_time/posix_time/posix_time.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        inline std::ostream& timestamp(std::ostream& s)
        {
            s << boost::posix_time::microsec_clock::local_time() << " ";
            return s;
        }
    }
}


#endif

