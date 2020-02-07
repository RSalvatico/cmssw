#ifndef CALOTSAMPLES_Pha2_H
#define CALOTSAMPLES_Pha2_H 1

#include "CalibFormats/CaloObjects/interface/CaloTSamplesBase.h"

/** \class CaloTSamples_Ph2

Class which represents the charge/voltage measurements of an event/channel
with the ADC decoding performed.

*/

template <class Ttype, uint32_t Tsize>
class CaloTSamples_Ph2 : public CaloTSamplesBase<Ttype> {
public:
  enum { kCapacity = Tsize };

  CaloTSamples_Ph2<Ttype, Tsize>();
  CaloTSamples_Ph2<Ttype, Tsize>(const CaloTSamples_Ph2<Ttype, Tsize> &cs);
  CaloTSamples_Ph2<Ttype, Tsize>(const DetId &id, uint32_t size = 0, uint32_t pre = 0);
  ~CaloTSamples_Ph2<Ttype, Tsize>() override;

  CaloTSamples_Ph2<Ttype, Tsize> &operator=(const CaloTSamples_Ph2<Ttype, Tsize> &cs);

  uint32_t capacity() const override;

private:
  Ttype *data(uint32_t i) override;
  const Ttype *cdata(uint32_t i) const override;

  Ttype m_data[Tsize];
};

#endif
