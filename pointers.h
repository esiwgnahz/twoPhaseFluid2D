#include <stdlib.h> // NULL

// base class for pointer and array classes
template<class T> class MCPtrBase
{
public:
  T&   operator  *() const { return *m_p; }
       operator T*() const { return m_p; }
  int  operator  !() const {return m_p==0;}
  T*   Relinquish(T* newp=0) { T* p = m_p; m_p = newp; return p; }

protected:
  //MCPtrBase(T* pointer) : m_p(pointer) {}
  //MCPtrBase() : m_p(0) {}

  T* m_p;

private:
  /*void* operator new(size_t);// {return 0;}  // prohibit use of new
  void operator delete[]( void * );
  void operator delete[]( void *, void * );
  //void operator delete( void * );
  void operator delete( void *, size_t );
  */
  /*
  #ifdef __GNUC__
  // in gcc operator delete seems to be needed for destructors in derived classes, so it can't be private
protected:
  #endif // __GNUC__
  void operator delete(void* p) {}// {((MCPtrBase<T>*)p)->m_p=0;}
  */
private:
  void operator++(); // declared but not implemented
  void operator++(int); // declared but not implemented
  void operator--(); // declared but not implemented
  void operator--(int); // declared but not implemented
};

// Pointer to a single object. Provides member access operator ->
// and uses object delete operator
template<class T> class MCPtr : public MCPtrBase<T>
{
public:
  /*explicit*/ MCPtr(T* pointer=NULL) { m_p=pointer; }
  ~MCPtr() { delete m_p; }
  MCPtr<T>& operator =(T* src) { delete m_p; m_p = src; return *this; }
  T* operator->() const { return m_p; }  // Could throw exception if m_p==0
private:
  #ifndef __BORLANDC__
  explicit MCPtr(const MCPtr<T> &src); // declared but not implemented
  #endif
  void operator =(const MCPtr<T> &src); // declared but not implemented
};

// Pointer to an array of type T. Provides an array subscript operator
// and uses array delete[]
template<class T> class MCAPtr : public MCPtrBase<T>
{
public:
  /*explicit*/ MCAPtr(T array[]=NULL) { m_p=array; }
  ~MCAPtr() { delete[] m_p; }
  MCAPtr<T>& operator =(T src[]) { delete[] m_p; m_p = src; return *this; }
  T& operator[](int i) const { return m_p[i]; }  // Could throw exception if m_p==0
  T& operator[](unsigned int i) const { return m_p[i]; }  // Could throw exception if m_p==0
  #ifdef _WIN64
  T& operator[](__int64 i) const { return m_p[i]; }  // Could throw exception if m_p==0
  T& operator[](unsigned __int64 i) const { return m_p[i]; }  // Could throw exception if m_p==0
  #endif
private: 
  #ifndef __BORLANDC__
  explicit MCAPtr(const MCAPtr<T> &src); // declared but not implemented
  #endif
  void operator =(const MCAPtr<T> &src); // declared but not implemented
};

template<class T> class MCRef : public MCPtrBase<T>
{
public:
  /*explicit*/ MCRef(T* pointer=NULL) { m_p=pointer; if(m_p) m_p->AddRef(); }
  ~MCRef() { if(m_p) m_p->Release(); }
  MCRef<T>& operator =(T* src) { if(m_p) m_p->Release(); m_p = src; if(m_p) m_p->AddRef(); return *this; }
  T* operator->() const { return m_p; }  // Could throw exception if m_p==0
  
  T** operator&() { if(m_p) { m_p->Release(); m_p=NULL; } return &m_p; }   //?? special operator! Use with care!
private:
  #ifndef __BORLANDC__
  explicit MCRef(const MCRef<T> &src); // declared but not implemented
  #endif
  void operator =(const MCRef<T> &src); // declared but not implemented
};

#ifndef _DEBUG
template<class T> void free(MCPtr<T>); // declared but not implemented
template<class T> void my_free(MCPtr<T>); // declared but not implemented
template<class T> void free(MCAPtr<T>); // declared but not implemented
template<class T> void my_free(MCAPtr<T>); // declared but not implemented
template<class T> void free(MCRef<T>); // declared but not implemented
template<class T> void my_free(MCRef<T>); // declared but not implemented
#endif


