#include "shastaLapack.hpp"

#include "algorithm.hpp"
#include "iostream.hpp"
#include "string.hpp"
#include "vector.hpp"

void shasta::testLapack()
{
    // Write out the copyright string so it is contained in the executable.
    cout << R"delimiter(
Shasta uses some functionality provided by the Lapack and Blas libraries.

Lapack and Blas come with this license:

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:
 .
 - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 .
 - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer listed
   in this license in the documentation and/or other materials
   provided with the distribution.
 .
 - Neither the name of the copyright holders nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
 .
 The copyright holders provide no reassurances that the source code
 provided does not infringe any patent, copyright, or any other
 intellectual property rights of third parties.  The copyright holders
 disclaim any liability to any recipient for claims brought against
 recipient by any third party for infringement of that parties
 intellectual property rights.
 .
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
)delimiter";


    // Solve the example here
    // https://en.wikipedia.org/wiki/Singular_value_decomposition
    // Compute a SVD of the following matrix:
    const string JOBU = "A";
    const string JOBVT = "A";
    const int M = 4;
    const int N = 5;
    vector<double> A = {
        1., 0., 0., 0.,
        0., 0., 0., 2.,
        0., 3., 0., 0.,
        0., 0., 0., 0.,
        2., 0., 0., 0.
    };
    const int LDA = M;
    vector<double> S(min(M, N));
    vector<double> U(M*M);
    const int LDU = M;
    vector<double> VT(N*N);
    const int LDVT = N;
    const int LWORK = 10 * max(M, N);
    vector<double> WORK(LWORK);
    int INFO = 0;

    dgesvd_(
        JOBU.data(), JOBVT.data(),
        M, N,
        &A[0], LDA, &S[0], &U[0], LDU, &VT[0], LDVT, &WORK[0], LWORK, INFO);
    cout << "dgesvd return code " << INFO << endl;

    if(INFO != 0) {
        return;
    }

    cout << "Singular values: " << endl;
    for(const double v: S) {
        cout << v << endl;
    }

}

