////////////////////////////////////////////////////////////////////////////
//
//  SortIndices.h
//      Defines a class that allows us to simultaneously sort C++ arrays.
//
////////////////////////////////////////////////////////////////////////////

class sort_indices
{
  private: 
    double * dparr;
  public:
    sort_indices(double*parr): dparr(parr){}
    bool operator()(int i, int j){ return dparr[i]>dparr[j]; }
};
