#ifndef ARRAY_INCLUDED 
#define ARRAY_INCLUDED

// copied shamelessly from Lippman
// chrisa 18-Nov-94


#ifdef SCCS_ID
static const char* arraySccsId = "@(#)array.h	1.4 09/15/95 11:43:21";
#endif


#include <iostream>

#include "assert.h"

const int iDefArraySize = 100; //10;  // if no arg to ctor this is size
const float fDefGrowthFactor = 1.5f;  // new size <-- old size multiplied by this
const float fAllocSlack = 1.05f;  // requested alloc <-- (needed * this_number)

template <class T>
class Array {
public:
  // ctor with no args gets empty array.  default ctor
  // for Arrays that are components of other classes.
  Array() : arrPtr(0), iSize(0), iUsed(0) {}

  // ctor with one int arg allocates but does not initialize
  // elements - only default constructor for T has been called.
  Array(const int iS) {
    int iNewSize = (int )(((float )iS * fAllocSlack) + 0.5);
    arrPtr = new T[iSize = iNewSize];
    iUsed = iS;  // though not initialized, they are available
  }

  // ctor can initialize from array of T and int size
  Array(const T *tyArr, const int iSize) : 
    arrPtr(0), iSize(0), iUsed(0) { init(tyArr, iSize); }

  // set array to a particular size.  deletes prev contents!
  void clearAndSetSize(const int iS) {
    // prev contents lost
    delete [] arrPtr;
    int iNewSize = (int )(((float )iS * fAllocSlack) + 0.5);
    arrPtr = new T[iSize = iNewSize];
    iUsed = iS;  // allow refs to these (empty) elems
  }
    
  // copy constructor.  size may change, since only iUsed are copied.
  Array(const Array<T> &arrRef) :
    arrPtr(0), iSize(0), iUsed(0) { init(arrRef.arrPtr, arrRef.iUsed); }

  ~Array() { delete [] arrPtr; }

  Array<T>& operator=(const Array<T>& aRef) {
    // never copy to yourself
    if (this != &aRef) {
      delete [] arrPtr;      // get rid of the old storage
      // re-allocate and copy the T refs stuff into yours
      // note this ignores actual storage size of passed ref
      init(aRef.arrPtr, aRef.iUsed);
    }
    return *this;      // and return ref to yourself
  }

  // returns the number of used elements, not storage size
  int length() const { return iUsed; }

  // insert a new element of T at pos
  // will bomb if you attempt to insert past first free
  // (i.e. iPos > length() )
  void insert(const int iPos, const T& newElement) {
    assert (iPos <= iUsed);
    if ((iUsed+1) > iSize) grow();  // no room? grow by default factor
    // shift contents of array from iPos to end down 1
    for (int i = iUsed-1; i >= iPos; i--) arrPtr[i+1] = arrPtr[i];
    // assign elem at iPos to passed ref
    arrPtr[iPos] = newElement;
    iUsed++;  // inc current size count
  }

  // insert at the end of the array
  void replace(const int iPos, const T& newElement) {
    assert (iPos <= iUsed);
    arrPtr[iPos] = newElement;
  }

  // insert at the end of the array
  void insertAtEnd(const T& newElement) {
    insert(length(), newElement);
  }

  // remove an element at pos, shift rest of array 
  void remove(const int iPos) {
    assert ((iPos >= 0) && (iPos < iUsed));
    for (int i = iPos; i < (iUsed-1); i++)
      arrPtr[i] = arrPtr[i+1];
    iUsed--;  // dec current size count, storage unchanged
  }

  // remove ALL elements that have the passed value
  // presumes == operator if T is a class
  void removeByValue(const T& killMe) {
    // read backward through the array
    for (int i = (iUsed - 1); i >= 0; i--) {
      if (arrPtr[i] == killMe) {
        remove(i);
      }
    }
  }

  // finds the FIRST element in the array whose value is
  // equal to that passed.  returns index if found, -1 if not
  int findByValue(const T& findMe) {
    for (int i = 0; i < iUsed; i++) {
      if (arrPtr[i] == findMe) { return(i); }
    }
    return (-1);
  }

  // an Array can be indexed into like a regular array
  T& operator [] (int iIndex) { 
    assert ((iIndex >= 0) && (iIndex < iUsed));
    return arrPtr[iIndex];
  }

  // same as above only returns const ref
  const T& operator [] (int iIndex) const { 
    assert ((iIndex >= 0) && (iIndex < iUsed));
    return arrPtr[iIndex];
  }

  // reverses the order of elements in the array, i.e.
  // a[0] <swapped with> a[iUsed], etc.
  void reverseOrder() {
    if (iUsed < 2) return;
    for (int i = 0; i < (iUsed / 2); i++) {
      T tTemp(arrPtr[i]);
      arrPtr[i] = arrPtr[iUsed - i - 1];
      arrPtr[iUsed -i -1] = tTemp;
    }
  }

private:
  // init the array of T with passed array
  void init (const T* array, const int iU) {
    // storage size is recalculated
    int iNewSize = (int )(((float )iU * fAllocSlack) + 0.5);
    arrPtr = new T[iSize = iNewSize];
    iUsed = iU;  // let them see what they passed
    assert ( arrPtr != 0 );
    for (int i = 0; i < iSize; i++) arrPtr[i] = array[i];
  }

  // allocs new mem of (factor * current) size, 
  // copies old to new.
  void grow(float fGrowthFactor = fDefGrowthFactor) {
    int iOldSize = iSize;   // save the old array size
    T *oldArrPtr = arrPtr;  // save the old array pointer
    int iNewSize;
    // increase the storage capacity by iGrowFactor
    if (iSize == 0) {
      iNewSize = iDefArraySize; // initially null array
    }
    else {
      iNewSize = (int )(((float )iOldSize * fGrowthFactor) + 0.5); 
    }
    arrPtr = new T[iSize = iNewSize];  assert(arrPtr != 0);
    // copy old to new
    for (int i = 0; i < iOldSize; i++) arrPtr[i] = oldArrPtr[i];
    delete [] oldArrPtr;    // delete the old storage
  }

  int iSize;        // size of currently allocated array
  int iUsed;        // how much is used (iUsed-1 == <index last elem>)
  T *arrPtr;     // ptr to array of T 
};

void foo();

#endif
   
