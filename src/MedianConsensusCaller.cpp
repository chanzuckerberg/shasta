#include "MedianConsensusCaller.hpp"
#include "Coverage.hpp"
#include <cmath>
using namespace ::shasta;


size_t MedianConsensusCaller::predict_runlength(const Coverage &coverage, AlignedBase consensus_base) const{
    size_t max_observed_repeat;     // Used to define the range of the loop over [base,length] observations
    size_t n_coverage;              // How many reads support each run length
    size_t n_total_coverage;        // How many reads support consensus base in total
    size_t prev_length;             // Which length was last visited in sorted observations
    size_t sum;                     // For tracking the cumulative coverage
    size_t median;                  // The median repeat length (this is returned)
    double midpoint;                // The point at which the CDF is 0.5

    max_observed_repeat = coverage.repeatCountEnd(consensus_base);
    n_total_coverage = coverage.coverage(consensus_base);

    midpoint = double(n_total_coverage)/2;

    sum = 0;
    median = 0;
    prev_length = 0;

    // Iterate possible observed lengths (including placeholders which may be 0 coverage)
    for (size_t length=0; length<=max_observed_repeat; length++){
        n_coverage = coverage.coverage(consensus_base, length);
        sum += n_coverage;

        if (double(sum) > midpoint){
            if (n_coverage > 1){
                // Both flanking observations are the same value
                median = length;
            }else{
                // Flanking observations differ in value (use the ceiling of their average)
                median = size_t(ceil(double(prev_length+length)/2));
            }
            break;
        }

        // Update temp variable for boundary cases where midpoint happens between 2 different observations
        if (n_coverage > 0){
            prev_length = length;
        }
    }

    return median;
}


Consensus MedianConsensusCaller::operator()(
    const Coverage& coverage) const
{
    const AlignedBase base = coverage.mostFrequentBase();
    const size_t repeatCount = predict_runlength(coverage, base);
    return Consensus(base, repeatCount);
}


void testMedianConsensusCaller(){
    MedianConsensusCaller classifier;
    Coverage coverage;
    AlignedBase consensus_base;
    Consensus consensus;

    // TEST CASE 1:
    // True index = 0.5
    // True median = 1
    // ignoring non-consensus bases

    // Arguments are base, strand, repeat count.
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 1);
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 1);
    coverage.addRead(AlignedBase::fromInteger((uint8_t)1), 0, 2);
    coverage.addRead(AlignedBase::fromInteger((uint8_t)2), 0, 2);
    coverage.addRead(AlignedBase::fromInteger((uint8_t)3), 0, 2);
    coverage.addRead(AlignedBase::fromInteger((uint8_t)4), 0, 0);

    consensus_base = coverage.mostFrequentBase();
    consensus = classifier(coverage);

    cout << "CONSENSUS BASE = " << consensus_base << "\n";
    cout << consensus.base << " " << consensus.repeatCount << "\n\n";

    // TEST CASE 2:
    // True index = 2.5
    // True median = 1

    coverage = Coverage();
    consensus_base = AlignedBase();

    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 1); // 1
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 1); // 2
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 1); // 3
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 1); // 4
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 2); // 5
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 3); // 6

    consensus_base = coverage.mostFrequentBase();
    consensus = classifier(coverage);

    cout << "CONSENSUS BASE = " << consensus_base << "\n";
    cout << consensus.base << " " << consensus.repeatCount << "\n\n";

    // TEST CASE 3:
    // True index = 3
    // True median = 2

    coverage = Coverage();
    consensus_base = AlignedBase();

    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 1); // 1
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 1); // 2
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 1); // 3
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 2); // 4
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 3); // 5
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 4); // 6
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 5); // 7

    consensus_base = coverage.mostFrequentBase();
    consensus = classifier(coverage);

    cout << "CONSENSUS BASE = " << consensus_base << "\n";
    cout << consensus.base << " " << consensus.repeatCount << "\n\n";

    // TEST CASE 4:
    // True index = 2.5
    // True median = 2

    coverage = Coverage();
    consensus_base = AlignedBase();

    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 1); // 1
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 1); // 2
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 1); // 3
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 3); // 4
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 4); // 5
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 5); // 6

    consensus_base = coverage.mostFrequentBase();
    consensus = classifier(coverage);

    cout << "CONSENSUS BASE = " << consensus_base << "\n";
    cout << consensus.base << " " << consensus.repeatCount << "\n\n";

    // TEST CASE 5
    // True index = 2.5
    // True median = 1.5
    // Rounded median = 2 (assuming ceiling rule)

    coverage = Coverage();
    consensus_base = AlignedBase();

    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 1); // 1
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 1); // 2
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 1); // 3
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 2); // 4
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 0, 3); // 5
    coverage.addRead(AlignedBase::fromInteger((uint8_t)0), 1, 4); // 6

    consensus_base = coverage.mostFrequentBase();
    consensus = classifier(coverage);

    cout << "CONSENSUS BASE = " << consensus_base << "\n";
    cout << consensus.base << " " << consensus.repeatCount << "\n\n";

}
