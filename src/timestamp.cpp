#include "timestamp.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "iostream.hpp"

std::ostream& shasta::timestamp(std::ostream& s)
{
    s << boost::posix_time::microsec_clock::local_time() << " ";
    return s;
}
