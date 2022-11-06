// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME WaveDataDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "inc/./WaveData.h"
#include "inc/./RunHeader.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_RunHeader(void *p = 0);
   static void *newArray_RunHeader(Long_t size, void *p);
   static void delete_RunHeader(void *p);
   static void deleteArray_RunHeader(void *p);
   static void destruct_RunHeader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RunHeader*)
   {
      ::RunHeader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RunHeader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RunHeader", ::RunHeader::Class_Version(), "inc/./RunHeader.h", 15,
                  typeid(::RunHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RunHeader::Dictionary, isa_proxy, 4,
                  sizeof(::RunHeader) );
      instance.SetNew(&new_RunHeader);
      instance.SetNewArray(&newArray_RunHeader);
      instance.SetDelete(&delete_RunHeader);
      instance.SetDeleteArray(&deleteArray_RunHeader);
      instance.SetDestructor(&destruct_RunHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RunHeader*)
   {
      return GenerateInitInstanceLocal((::RunHeader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RunHeader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WaveData(void *p = 0);
   static void *newArray_WaveData(Long_t size, void *p);
   static void delete_WaveData(void *p);
   static void deleteArray_WaveData(void *p);
   static void destruct_WaveData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WaveData*)
   {
      ::WaveData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WaveData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WaveData", ::WaveData::Class_Version(), "inc/./WaveData.h", 12,
                  typeid(::WaveData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WaveData::Dictionary, isa_proxy, 4,
                  sizeof(::WaveData) );
      instance.SetNew(&new_WaveData);
      instance.SetNewArray(&newArray_WaveData);
      instance.SetDelete(&delete_WaveData);
      instance.SetDeleteArray(&deleteArray_WaveData);
      instance.SetDestructor(&destruct_WaveData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WaveData*)
   {
      return GenerateInitInstanceLocal((::WaveData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::WaveData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RunHeader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RunHeader::Class_Name()
{
   return "RunHeader";
}

//______________________________________________________________________________
const char *RunHeader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RunHeader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RunHeader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RunHeader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RunHeader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RunHeader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RunHeader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RunHeader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WaveData::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WaveData::Class_Name()
{
   return "WaveData";
}

//______________________________________________________________________________
const char *WaveData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WaveData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WaveData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WaveData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WaveData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WaveData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WaveData::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WaveData*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RunHeader::Streamer(TBuffer &R__b)
{
   // Stream an object of class RunHeader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RunHeader::Class(),this);
   } else {
      R__b.WriteClassBuffer(RunHeader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RunHeader(void *p) {
      return  p ? new(p) ::RunHeader : new ::RunHeader;
   }
   static void *newArray_RunHeader(Long_t nElements, void *p) {
      return p ? new(p) ::RunHeader[nElements] : new ::RunHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_RunHeader(void *p) {
      delete ((::RunHeader*)p);
   }
   static void deleteArray_RunHeader(void *p) {
      delete [] ((::RunHeader*)p);
   }
   static void destruct_RunHeader(void *p) {
      typedef ::RunHeader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RunHeader

//______________________________________________________________________________
void WaveData::Streamer(TBuffer &R__b)
{
   // Stream an object of class WaveData.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WaveData::Class(),this);
   } else {
      R__b.WriteClassBuffer(WaveData::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WaveData(void *p) {
      return  p ? new(p) ::WaveData : new ::WaveData;
   }
   static void *newArray_WaveData(Long_t nElements, void *p) {
      return p ? new(p) ::WaveData[nElements] : new ::WaveData[nElements];
   }
   // Wrapper around operator delete
   static void delete_WaveData(void *p) {
      delete ((::WaveData*)p);
   }
   static void deleteArray_WaveData(void *p) {
      delete [] ((::WaveData*)p);
   }
   static void destruct_WaveData(void *p) {
      typedef ::WaveData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WaveData

namespace ROOT {
   static TClass *vectorlEunsignedsPshortgR_Dictionary();
   static void vectorlEunsignedsPshortgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPshortgR(void *p = 0);
   static void *newArray_vectorlEunsignedsPshortgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPshortgR(void *p);
   static void deleteArray_vectorlEunsignedsPshortgR(void *p);
   static void destruct_vectorlEunsignedsPshortgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned short>*)
   {
      vector<unsigned short> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned short>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned short>", -2, "vector", 216,
                  typeid(vector<unsigned short>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPshortgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<unsigned short>) );
      instance.SetNew(&new_vectorlEunsignedsPshortgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPshortgR);
      instance.SetDelete(&delete_vectorlEunsignedsPshortgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPshortgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPshortgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned short> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<unsigned short>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPshortgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<unsigned short>*)0x0)->GetClass();
      vectorlEunsignedsPshortgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPshortgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPshortgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned short> : new vector<unsigned short>;
   }
   static void *newArray_vectorlEunsignedsPshortgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned short>[nElements] : new vector<unsigned short>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPshortgR(void *p) {
      delete ((vector<unsigned short>*)p);
   }
   static void deleteArray_vectorlEunsignedsPshortgR(void *p) {
      delete [] ((vector<unsigned short>*)p);
   }
   static void destruct_vectorlEunsignedsPshortgR(void *p) {
      typedef vector<unsigned short> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned short>

namespace ROOT {
   static TClass *vectorlEunsignedsPintgR_Dictionary();
   static void vectorlEunsignedsPintgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPintgR(void *p = 0);
   static void *newArray_vectorlEunsignedsPintgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPintgR(void *p);
   static void deleteArray_vectorlEunsignedsPintgR(void *p);
   static void destruct_vectorlEunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned int>*)
   {
      vector<unsigned int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned int>", -2, "vector", 216,
                  typeid(vector<unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPintgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<unsigned int>) );
      instance.SetNew(&new_vectorlEunsignedsPintgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPintgR);
      instance.SetDelete(&delete_vectorlEunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPintgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<unsigned int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<unsigned int>*)0x0)->GetClass();
      vectorlEunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int> : new vector<unsigned int>;
   }
   static void *newArray_vectorlEunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int>[nElements] : new vector<unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPintgR(void *p) {
      delete ((vector<unsigned int>*)p);
   }
   static void deleteArray_vectorlEunsignedsPintgR(void *p) {
      delete [] ((vector<unsigned int>*)p);
   }
   static void destruct_vectorlEunsignedsPintgR(void *p) {
      typedef vector<unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned int>

namespace ROOT {
   static TClass *unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR_Dictionary();
   static void unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR_TClassManip(TClass*);
   static void *new_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR(void *p = 0);
   static void *newArray_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR(Long_t size, void *p);
   static void delete_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR(void *p);
   static void deleteArray_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR(void *p);
   static void destruct_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<unsigned short,vector<unsigned short> >*)
   {
      unordered_map<unsigned short,vector<unsigned short> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<unsigned short,vector<unsigned short> >));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<unsigned short,vector<unsigned short> >", -2, "unordered_map", 101,
                  typeid(unordered_map<unsigned short,vector<unsigned short> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(unordered_map<unsigned short,vector<unsigned short> >) );
      instance.SetNew(&new_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR);
      instance.SetNewArray(&newArray_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR);
      instance.SetDelete(&delete_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR);
      instance.SetDestructor(&destruct_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<unsigned short,vector<unsigned short> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const unordered_map<unsigned short,vector<unsigned short> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const unordered_map<unsigned short,vector<unsigned short> >*)0x0)->GetClass();
      unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) unordered_map<unsigned short,vector<unsigned short> > : new unordered_map<unsigned short,vector<unsigned short> >;
   }
   static void *newArray_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) unordered_map<unsigned short,vector<unsigned short> >[nElements] : new unordered_map<unsigned short,vector<unsigned short> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR(void *p) {
      delete ((unordered_map<unsigned short,vector<unsigned short> >*)p);
   }
   static void deleteArray_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR(void *p) {
      delete [] ((unordered_map<unsigned short,vector<unsigned short> >*)p);
   }
   static void destruct_unordered_maplEunsignedsPshortcOvectorlEunsignedsPshortgRsPgR(void *p) {
      typedef unordered_map<unsigned short,vector<unsigned short> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class unordered_map<unsigned short,vector<unsigned short> >

namespace {
  void TriggerDictionaryInitialization_WaveDataDict_Impl() {
    static const char* headers[] = {
"inc/./WaveData.h",
"inc/./RunHeader.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf7/include",
"/lustrehome/bdanzi/offline_analysis/dchdatareade_test_beam/WaveData/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "WaveDataDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$inc/./WaveData.h")))  RunHeader;
class __attribute__((annotate("$clingAutoload$inc/./WaveData.h")))  WaveData;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "WaveDataDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "inc/./WaveData.h"
#include "inc/./RunHeader.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RunHeader", payloadCode, "@",
"WaveData", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("WaveDataDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_WaveDataDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_WaveDataDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_WaveDataDict() {
  TriggerDictionaryInitialization_WaveDataDict_Impl();
}
