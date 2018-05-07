#include "LongBaseSequence.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

#include "vector.hpp"



void LongBaseSequences::createNew(
    const string& name,
    size_t pageSize)
{
    baseCount.createNew(name + "-BaseCount", pageSize);
    data.createNew(name + "-Bases", pageSize);
}



void LongBaseSequences::accessExistingReadOnly(const string& name)
{
    baseCount.accessExistingReadOnly(name + "-BaseCount");
    data.accessExistingReadOnly(name + "-Bases");
    CZI_ASSERT(baseCount.size() == data.size());
}



void LongBaseSequences::accessExistingReadWrite(const string& name)
{
    baseCount.accessExistingReadWrite(name + "-BaseCount");
    data.accessExistingReadWrite(name + "-Bases");
    CZI_ASSERT(baseCount.size() == data.size());
}



void LongBaseSequences::accessExistingReadWriteOrCreateNew(
    const string& name,
    size_t pageSize)
{
    baseCount.accessExistingReadWriteOrCreateNew(name + "-BaseCount", pageSize);
    data.accessExistingReadWriteOrCreateNew(name + "-Bases", pageSize);
    CZI_ASSERT(baseCount.size() == data.size());
}



void LongBaseSequences::remove()
{
    baseCount.remove();
    data.remove();
}


// Append a new sequence at the end.
void LongBaseSequences::append(const LongBaseSequence& s)
{
    baseCount.push_back(s.baseCount);
    data.appendVector(s.data.begin(), s.data.end());
}
void LongBaseSequences::append(const vector<Base>& s)
{
    append(LongBaseSequence(s));
}
void LongBaseSequences::append(size_t baseCountArgument)
{
    baseCount.push_back(baseCountArgument);
    const size_t wordCount = LongBaseSequenceView::wordCount(baseCountArgument);
    data.appendVector(wordCount);
    CZI_ASSERT(baseCount.size() == data.size());
}



void ChanZuckerberg::Nanopore2::testLongBaseSequence()
{



    // Test class LongBaseSequenceView.
    {
        vector<uint64_t> v(4, 0);
        LongBaseSequenceView s(v.data(), 100);
        s.set( 5, Base('C', Base::FromCharacter()));
        s.set(63, Base('G', Base::FromCharacter()));
        s.set(64, Base('T', Base::FromCharacter()));
        s.set(95, Base('T', Base::FromCharacter()));

        for(size_t i=0; i<100; i++) {
            cout << s[i];
        }
        cout << endl;

        /*
        const auto oldFill = cout.fill('0');
        for(const uint64_t x: v) {
            cout << std::setw(16) << std::hex << x << endl;
        }
        cout.fill(oldFill);
        */
    }



    // Test class LongBaseSequence.
    {
        LongBaseSequence s(100);
        s.set( 5, Base('C', Base::FromCharacter()));
        s.set(63, Base('G', Base::FromCharacter()));
        s.set(64, Base('T', Base::FromCharacter()));
        s.set(95, Base('T', Base::FromCharacter()));

        for(size_t i=0; i<100; i++) {
            cout << s[i];
        }
        cout << endl;

        /*
        const auto oldFill = cout.fill('0');
        for(const uint64_t x: v) {
            cout << sad::setw(16) << std::hex << x << endl;
        }
        cout.fill(oldFill);
        */
    }



    // Test class LongBaseSequences.
    {
        LongBaseSequences sequences;
        sequences.createNew("abc", 4096);

        for(size_t i=0; i<8; i++) {
            LongBaseSequence sequence(i*50);
            const Base base(i%4, Base::FromInteger());
            for(size_t j=0; j<sequence.baseCount; j++) {
                sequence.set(j, base);
            }
            cout << sequence << endl;
            sequences.append(sequence);
        }
        cout << endl;
        for(size_t i=0; i<sequences.size(); i++) {\
            cout << i << endl;
            LongBaseSequenceView sequence = sequences[i];
            const Base base(i%4, Base::FromInteger());
            for(size_t j=0; j<sequence.baseCount; j++) {
                sequence.set(j, base);
            }
            cout << sequence << endl;
        }
    }
}
