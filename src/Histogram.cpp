#include "Histogram.hpp"
#include <unordered_set>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace shasta;
using std::unordered_set;
using std::runtime_error;
using std::to_string;
using std::string;

shasta::Histogram2::Histogram2(
        double start,
        double stop,
        size_t binCount,
        bool unboundedLeft,
        bool unboundedRight):
        start(start),
        stop(stop),
        binCount(binCount),
        binSize((stop-start)/double(binCount)),
        histogram(binCount, 0),
        unboundedLeft(unboundedLeft),
        unboundedRight(unboundedRight)
{}


int64_t shasta::Histogram2::findIndex(double x){
    /// Find index of bin by normalizing the value w.r.t. bin edges.
    /// If the value does not fit in the range of the histogram, this function returns -1
    /// If the user wants to accumulate values outside the bounds defined, then any value out of bounds will result in
    /// the boundary bins being incremented
    auto index = int64_t(floor((x - start) / binSize));

    if (x < start){
        if(unboundedLeft) {
            index = 0;
        }
        else{
            index = -1;
        }
    }

    if (x >= stop) {
        if (unboundedRight) {
            index = int64_t(binCount) - 1;
        }
        else{
            index = -1;
        }
    }

    return index;
}


void shasta::Histogram2::update(double x) {
    auto index = findIndex(x);

    if (index >= 0) {
        histogram[size_t(index)]++;
    }
    else {
        // Do nothing. Print warning?
//        std::cerr << "Value " << x << " ignored because index is " << index << ". Size: " << histogram.size() << '\n';
    }
}


void shasta::Histogram2::getNormalizedHistogram(vector<double>& normalizedHistogram){
    uint64_t sum = getSum();

    for (auto& e: histogram){
        normalizedHistogram.push_back(double(e)/double(sum));
    }
}


uint64_t Histogram2::getSum(){
    uint64_t sum = 0;
    for (auto& item: histogram){
        sum += item;
    }
    return sum;
}


void shasta::Histogram2::writeToHtml(ostream& html, uint64_t sizePx){
    uint64_t yMax = 0;
    for (auto& e: histogram){
        if (e > yMax){
            yMax = e;
        }
    }

    double scale = double(sizePx)/double(yMax);

    html << "<table style='margin-top: 1em; margin-bottom: 1em'>";
    html << "<tr>"
            "<th class='centered'>Left bound"
            "<th class='centered'>Right bound"
            "<th class='centered'>Count"
            "<th class='centered'>Plot";

    int32_t precision;

    for (size_t i=0; i<histogram.size(); i++){
        const double leftBound = double(i)*(binSize);
        const double rightBound = double(i+1)*(binSize);
        const auto y = histogram[i];

        // Check if the bin bounds have any trailing decimals
        if (std::fmod(leftBound,1) == 0 and std::fmod(rightBound,1) == 0){
            precision = 0;
        }
        else{
            precision = 2;
        }

        string leftBoundString;
        string rightBoundString;

        if (unboundedLeft and i==0){
            leftBoundString = "-inf";
        }
        else {
            leftBoundString = to_string(leftBound);
            const size_t decimalPosition = leftBoundString.find('.');
            leftBoundString = leftBoundString.substr(0,decimalPosition+3);
        }

        if (unboundedRight and i == binCount-1){
            rightBoundString = "inf";
        }
        else{
            rightBoundString = to_string(rightBound);
            const size_t decimalPosition = rightBoundString.find('.');
            rightBoundString = rightBoundString.substr(0,decimalPosition+3);

        }


        html << std::fixed << std::setprecision(precision) <<
             "<tr>"
             "<td class=centered>" << leftBoundString <<
             "<td class=centered>" << rightBoundString <<
             "<td class=centered>" << y <<
             "<td>"
             "<div class=sketch title='alignedFractionHistogram' style='display:inline-block;margin:0px;padding:0px;"
             "background-color:blue;height:6px;width:" << double(y)*scale << "px;'></div>";
    }
    html << "</table>";

    // Remove precision settings that were specified above
    html.unsetf(std::ios_base::floatfield);
}


void shasta::writeHistogramsToHtml(ostream& html, Histogram2& histogramA, Histogram2& histogramB, uint64_t sizePx){
    // First verify that histograms are compatible
    bool compatible = true;
    if (histogramA.start != histogramB.start){
        compatible = false;
    }
    if (histogramA.stop != histogramB.stop){
        compatible = false;
    }
    if (histogramA.binCount != histogramB.binCount){
        compatible = false;
    }
    if (histogramA.unboundedLeft != histogramB.unboundedLeft){
        compatible = false;
    }
    if (histogramA.unboundedRight != histogramB.unboundedRight){
        compatible = false;
    }

    if (not compatible){
        throw runtime_error("ERROR: histograms with differing bins cannot be plotted together");
    }

    // Find max frequency, for plot scaling purposes
    uint64_t yMax = 0;
    for (size_t i=0; i<histogramA.histogram.size(); i++){
        if (histogramA.histogram[i] > yMax) {
            yMax = histogramA.histogram[i];
        }
        if (histogramB.histogram[i] > yMax) {
            yMax = histogramB.histogram[i];
        }
    }

    double scale = double(sizePx)/double(yMax);

    html << "<table style='margin-top: 1em; margin-bottom: 1em'>";
    html << "<tr>"
            "<th class='centered'>Left bound"
            "<th class='centered'>Right bound"
            "<th class='centered'>Count A"
            "<th class='centered'>Count B"
            "<th class='centered'>Plot";

    int32_t precision;

    for (size_t i=0; i<histogramA.histogram.size(); i++){
        const double leftBound = double(i)*(histogramA.binSize);
        const double rightBound = double(i+1)*(histogramA.binSize);

        // Check if the bin bounds have any trailing decimals
        if (std::fmod(leftBound,1) == 0 and std::fmod(rightBound,1) == 0){
            precision = 0;
        }
        else{
            precision = 2;
        }

        string leftBoundString;
        string rightBoundString;

        if (histogramA.unboundedLeft and i==0){
            leftBoundString = "-inf";
        }
        else {
            leftBoundString = to_string(leftBound);
            const size_t decimalPosition = leftBoundString.find('.');
            leftBoundString = leftBoundString.substr(0,decimalPosition+3);
        }

        if (histogramA.unboundedRight and i == histogramA.binCount-1){
            rightBoundString = "inf";
        }
        else{
            rightBoundString = to_string(rightBound);
            const size_t decimalPosition = rightBoundString.find('.');
            rightBoundString = rightBoundString.substr(0,decimalPosition+3);

        }

        html << std::fixed << std::setprecision(precision) <<
            "<tr>"
            "<td class=centered>" << leftBoundString <<
            "<td class=centered>" << rightBoundString <<
            "<td class=centered>" << histogramA.histogram[i] <<
            "<td class=centered>" << histogramB.histogram[i] <<
            "<td style='line-height:8px;white-space:nowrap'>" <<
            "<div class=sketch title='alignedFractionHistogram' "
            "style='display:inline-block;margin:0px;padding:0px;background-color:red;height:6px;width:" <<
            double(histogramA.histogram[i]) * scale << "px;'></div><br>" <<
            "<div class=sketch title='alignedFractionHistogram' "
            "style='display:inline-block;margin:0px;padding:0px;background-color:blue;height:6px;width:" <<
            double(histogramB.histogram[i]) * scale << "px;'></div>";
    }
    html << "</table>";

    // Remove precision settings that were specified above
    html.unsetf(std::ios_base::floatfield);
}


void shasta::testIterativeHistogram(){
    /// Test the iterative histogram

    std::cout << "\nTESTING Histogram2 case 1\n";

    Histogram2 h0(0, 10, 10);
    h0.update(0);        // 1
    h0.update(-1);       // None
    h0.update(10);       // None
    h0.update(9.99999);  // 10
    h0.update(10.0001);  // None
    h0.update(0.5);      // 1
    h0.update(1.5);      // 2
    h0.update(1.0);      // 2
    h0.update(1.99999);  // 2

    // ^ expect [2,3,0,0,0,0,0,0,0,1]

    for (auto& e: h0.histogram) {
        std::cout << e << ",";
    }
    std::cout << '\n';

    std::cout << "\nTESTING Histogram2 case 2\n";

    Histogram2 h1(0, 1.0, 10);
    h1.update(0);         // 1
    h1.update(-0.1);      // None
    h1.update(1.0);       // None
    h1.update(0.999999);  // 10
    h1.update(1.00001);   // None
    h1.update(0.05);      // 1
    h1.update(0.15);      // 2
    h1.update(0.10);      // 2
    h1.update(0.199999);  // 2

    // ^ expect [2,3,0,0,0,0,0,0,0,1]

    for (auto& e: h1.histogram) {
        std::cout << e << ",";
    }
    std::cout << '\n';

    std::cout << "\nTESTING Histogram2 case 3\n";

    Histogram2 h2(1, 2.0, 10);
    h2.update(1 + 0);         // 1
    h2.update(1 + -0.1);      // None
    h2.update(1 + 1.0);       // None
    h2.update(1 + 0.999999);  // 10
    h2.update(1 + 1.00001);   // None
    h2.update(1 + 0.05);      // 1
    h2.update(1 + 0.15);      // 2
    h2.update(1 + 0.10);      // 2
    h2.update(1 + 0.199999);  // 2

    // ^ expect [2,3,0,0,0,0,0,0,0,1]

    for (auto& e: h2.histogram) {
        std::cout << e << ",";
    }
    std::cout << '\n';

    std::cout << "\nTESTING Histogram2 case 4\n";

    Histogram2 h3(-0.5, 0.5, 10);
    h3.update(-0.5 + 0);         // 1
    h3.update(-0.5 + -0.1);      // None
    h3.update(-0.5 + 1.0);       // None right edge
    h3.update(-0.5 + 0.999999);  // 10
    h3.update(-0.5 + 1.00001);   // None
    h3.update(-0.5 + 0.05);      // 1
    h3.update(-0.5 + 0.15);      // 2
    h3.update(-0.5 + 0.10);      // 2 ... in near-edge cases float division may shift left a bin
    h3.update(-0.5 + 0.199999);  // 2

    // ^ expect [2,3,0,0,0,0,0,0,0,1]

    for (auto& e: h3.histogram) {
        std::cout << e << ",";
    }
    std::cout << '\n';

    vector<double> normalized3;
    h3.getNormalizedHistogram(normalized3);
    double sum3 = 0;

    for (auto& e: normalized3) {
        std::cout << e << ",";
        sum3 += e;
    }
    std::cout << '\n';
    std::cout << sum3 << '\n';

    std::cout << "\nTESTING Histogram2 case 5\n";

    Histogram2 h4(0, 1.0, 10, true, true);
    h4.update(0);         // 1
    h4.update(-0.1);      // 1
    h4.update(1.0);       // 10
    h4.update(0.999999);  // 10
    h4.update(1.00001);   // 10
    h4.update(0.05);      // 1
    h4.update(0.15);      // 2
    h4.update(0.10);      // 2
    h4.update(0.199999);  // 2

    // ^ expect [3,3,0,0,0,0,0,0,0,3]

    for (auto& e: h4.histogram) {
        std::cout << e << ",";
    }
    std::cout << '\n';

    vector<double> normalized4;
    h4.getNormalizedHistogram(normalized4);
    double sum4 = 0;

    for (auto& e: normalized4) {
        std::cout << e << ",";
        sum4 += e;
    }
    std::cout << '\n';
    std::cout << sum4 << '\n';

}