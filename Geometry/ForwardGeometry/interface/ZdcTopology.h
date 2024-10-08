#ifndef Geometry_ForwardGeometry_ZdcTopology_H
#define Geometry_ForwardGeometry_ZdcTopology_H 1

#include <vector>
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/HcalCommonData/interface/HcalTopologyMode.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"

/** \class ZDCTopology

   $Revision: 1.0 $
   \author E. Garcia - UIC
*/

class ZdcTopology : public CaloSubdetectorTopology {
public:
  ZdcTopology(const HcalDDDRecConstants* hcons);

  HcalTopologyMode::Mode mode() const { return mode_; }
  /** Exlucde a cell*/
  void exclude(const HcalZDCDetId& id);
  /** Exclude a side*/
  void exclude(int zside);
  /** Exclude a section, in either side (+1 positive, -1 negative)*/
  void exclude(int zside, HcalZDCDetId::Section section);
  /** Exclude a range of channels (deph) for a given subdetector*/
  int exclude(int zside, HcalZDCDetId::Section section, int ich1, int ich2);

  /** Is this a valid cell id? */
  using CaloSubdetectorTopology::valid;
  virtual bool valid(const HcalZDCDetId& id) const;

  /** Get the transverse (X) neighbors of the given cell*/
  virtual std::vector<DetId> transverse(const DetId& id) const;

  /** Get the longitudinal neighbors (Z) of the given cell*/
  virtual std::vector<DetId> longitudinal(const DetId& id) const;

  //** I have to put this here since they inherit from CaloSubdetectorTopology
  std::vector<DetId> east(const DetId& id) const override;
  std::vector<DetId> west(const DetId& id) const override;
  std::vector<DetId> north(const DetId& id) const override;
  std::vector<DetId> south(const DetId& id) const override;
  std::vector<DetId> up(const DetId& id) const override;
  std::vector<DetId> down(const DetId& id) const override;

  // how many channels (deph) for a given section
  using CaloSubdetectorTopology::ncells;
  int ncells(HcalZDCDetId::Section section) const;

  //return first and last cell of each section
  int firstCell(HcalZDCDetId::Section section) const;
  int lastCell(HcalZDCDetId::Section section) const;

  uint32_t kSizeForDenseIndexing() const;
  bool validDenseIndex(uint32_t din) const { return (din < kSizeForDenseIndexing()); }

  DetId denseId2detId(uint32_t di) const override;
  uint32_t detId2DenseIndex(const DetId& id) const;

private:
  bool validRaw(const HcalZDCDetId& id) const;

  bool isExcluded(const HcalZDCDetId& id) const;

  int firstEMModule() const { return firstEMModule_; }
  int firstHADModule() const { return firstHADModule_; }
  int firstLUMModule() const { return firstLUMModule_; }
  int firstRPDModule() const { return firstRPDModule_; }
  int lastEMModule() const { return lastEMModule_; }
  int lastHADModule() const { return lastHADModule_; }
  int lastLUMModule() const { return lastLUMModule_; }
  int lastRPDModule() const { return lastRPDModule_; }

  const HcalDDDRecConstants* hcons_;
  HcalTopologyMode::Mode mode_;

  std::vector<HcalZDCDetId> exclusionList_;

  bool excludeEM_, excludeHAD_, excludeLUM_, excludeRPD_, excludeZP_, excludeZN_;

  int firstEMModule_, lastEMModule_, firstHADModule_, lastHADModule_, firstLUMModule_, lastLUMModule_, firstRPDModule_,
      lastRPDModule_;
};

#endif
