#include "LongBaseSequence.hpp"
using namespace shasta;

#include "vector.hpp"



void LongBaseSequences::createNew(
    const string& name,
    size_t pageSize)
{
    if(name.empty()) {
        baseCount.createNew("", pageSize);
        data.createNew("", pageSize);
    } else {
        baseCount.createNew(name + "-BaseCount", pageSize);
        data.createNew(name + "-Bases", pageSize);
    }
}



void LongBaseSequences::accessExistingReadOnly(const string& name)
{
    baseCount.accessExistingReadOnly(name + "-BaseCount");
    data.accessExistingReadOnly(name + "-Bases");
    SHASTA_ASSERT(baseCount.size() == data.size());
}



void LongBaseSequences::accessExistingReadWrite(const string& name)
{
    baseCount.accessExistingReadWrite(name + "-BaseCount");
    data.accessExistingReadWrite(name + "-Bases");
    SHASTA_ASSERT(baseCount.size() == data.size());
}



void LongBaseSequences::accessExistingReadWriteOrCreateNew(
    const string& name,
    size_t pageSize)
{
    baseCount.accessExistingReadWriteOrCreateNew(name + "-BaseCount", pageSize);
    data.accessExistingReadWriteOrCreateNew(name + "-Bases", pageSize);
    SHASTA_ASSERT(baseCount.size() == data.size());
}



void LongBaseSequences::remove()
{
    baseCount.remove();
    data.remove();
}
void LongBaseSequences::close()
{
    baseCount.close();
    data.close();
}

void LongBaseSequences::clear()
{
    baseCount.clear();
    data.clear();
}

void LongBaseSequences::unreserve() {
    baseCount.unreserve();
    data.unreserve();
}

// Append a new sequence at the end.
void LongBaseSequences::append(const LongBaseSequenceView& s)
{
    baseCount.push_back(s.baseCount);
    data.appendVector(s.begin, s.begin+LongBaseSequenceView::wordCount(s.baseCount));
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
    SHASTA_ASSERT(baseCount.size() == data.size());
}



void shasta::testLongBaseSequence()
{



    // Test class LongBaseSequenceView.
    {
        vector<uint64_t> v(4, 0);
        LongBaseSequenceView s(v.data(), 100);
        s.set( 5, Base::fromCharacter('C'));
        s.set(63, Base::fromCharacter('G'));
        s.set(64, Base::fromCharacter('T'));
        s.set(95, Base::fromCharacter('T'));

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
        s.set( 5, Base::fromCharacter('C'));
        s.set(63, Base::fromCharacter('G'));
        s.set(64, Base::fromCharacter('T'));
        s.set(95, Base::fromCharacter('T'));

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

        for(uint8_t i=0; i<8; i++) {
            LongBaseSequence sequence(i*50);
            const Base base = Base::fromInteger(uint8_t(i%4));
            for(size_t j=0; j<sequence.baseCount; j++) {
                sequence.set(j, base);
            }
            cout << sequence << endl;
            sequences.append(sequence);
        }
        cout << endl;
        for(uint8_t i=0; i<sequences.size(); i++) {\
            cout << i << endl;
            LongBaseSequenceView sequence = sequences[i];
            const Base base = Base::fromInteger(uint8_t(i%4));
            for(size_t j=0; j<sequence.baseCount; j++) {
                sequence.set(j, base);
            }
            cout << sequence << endl;
        }
    }
}
