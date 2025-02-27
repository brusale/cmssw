#ifndef FWCore_Framework_EventPrincipal_h
#define FWCore_Framework_EventPrincipal_h

/*----------------------------------------------------------------------

EventPrincipal: This is the class responsible for management of
per event EDProducts. It is not seen by reconstruction code;
such code sees the Event class, which is a proxy for EventPrincipal.

The major internal component of the EventPrincipal
is the DataBlock.

----------------------------------------------------------------------*/

#include "DataFormats/Common/interface/WrapperBase.h"
#include "DataFormats/Provenance/interface/BranchListIndex.h"
#include "FWCore/Framework/interface/ProductProvenanceRetriever.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Provenance/interface/EventSelectionID.h"
#include "DataFormats/Provenance/interface/EventToProcessBlockIndexes.h"
#include "FWCore/Common/interface/FWCoreCommonFwd.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/Signal.h"
#include "FWCore/Utilities/interface/get_underlying_safe.h"
#include "FWCore/Framework/interface/Principal.h"
#include "FWCore/Framework/interface/ProductResolversFactory.h"

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <optional>

namespace edm {
  class BranchID;
  class BranchIDListHelper;
  class ProductProvenanceRetriever;
  class DelayedReader;
  class EventID;
  class HistoryAppender;
  class LuminosityBlockPrincipal;
  class ModuleCallingContext;
  class ProductID;
  class StreamContext;
  class ThinnedAssociation;
  class ThinnedAssociationsHelper;
  class RunPrincipal;

  class EventPrincipal : public Principal {
  public:
    typedef EventAuxiliary Auxiliary;
    typedef Principal Base;

    typedef Base::ConstProductResolverPtr ConstProductResolverPtr;
    static constexpr int invalidBunchXing = EventAuxiliary::invalidBunchXing;
    static constexpr int invalidStoreNumber = EventAuxiliary::invalidStoreNumber;
    template <ProductResolversFactory FACTORY>
    EventPrincipal(std::shared_ptr<ProductRegistry const> reg,
                   FACTORY&& iFactory,
                   std::shared_ptr<BranchIDListHelper const> branchIDListHelper,
                   std::shared_ptr<ThinnedAssociationsHelper const> thinnedAssociationsHelper,
                   ProcessConfiguration const& pc,
                   HistoryAppender* historyAppender,
                   unsigned int streamIndex = 0,
                   ProcessBlockHelperBase const* processBlockHelper = nullptr)
        : EventPrincipal(reg,
                         iFactory(InEvent, pc.processName(), *reg),
                         branchIDListHelper,
                         thinnedAssociationsHelper,
                         pc,
                         historyAppender,
                         streamIndex,
                         processBlockHelper) {}

    ~EventPrincipal() override {}

    void fillEventPrincipal(EventAuxiliary const& aux,
                            ProcessHistory const* processHistory,
                            DelayedReader* reader = nullptr);
    void fillEventPrincipal(EventAuxiliary const& aux,
                            ProcessHistory const* processHistory,
                            EventSelectionIDVector eventSelectionIDs,
                            BranchListIndexes branchListIndexes,
                            DelayedReader* reader = nullptr);
    //provRetriever is changed via a call to ProductProvenanceRetriever::deepSwap
    void fillEventPrincipal(EventAuxiliary const& aux,
                            ProcessHistory const* processHistory,
                            EventSelectionIDVector eventSelectionIDs,
                            BranchListIndexes branchListIndexes,
                            EventToProcessBlockIndexes const&,
                            ProductProvenanceRetriever const& provRetriever,
                            DelayedReader* reader = nullptr,
                            bool deepCopyRetriever = true);

    void clearEventPrincipal();

    LuminosityBlockPrincipal const& luminosityBlockPrincipal() const { return *luminosityBlockPrincipal_; }

    LuminosityBlockPrincipal& luminosityBlockPrincipal() { return *luminosityBlockPrincipal_; }

    bool luminosityBlockPrincipalPtrValid() const { return luminosityBlockPrincipal_ != nullptr; }

    //does not share ownership
    void setLuminosityBlockPrincipal(LuminosityBlockPrincipal* lbp);

    void setRunAndLumiNumber(RunNumber_t run, LuminosityBlockNumber_t lumi);

    EventID const& id() const { return aux().id(); }

    Timestamp const& time() const { return aux().time(); }

    bool isReal() const { return aux().isRealData(); }

    EventAuxiliary::ExperimentType ExperimentType() const { return aux().experimentType(); }

    int bunchCrossing() const { return aux().bunchCrossing(); }

    int storeNumber() const { return aux().storeNumber(); }

    EventAuxiliary const& aux() const { return aux_; }

    StreamID streamID() const { return streamID_; }

    LuminosityBlockNumber_t luminosityBlock() const { return id().luminosityBlock(); }

    RunNumber_t run() const { return id().run(); }

    RunPrincipal const& runPrincipal() const;

    ProductProvenanceRetriever const* productProvenanceRetrieverPtr() const { return provRetrieverPtr_.get(); }

    EventSelectionIDVector const& eventSelectionIDs() const;

    BranchListIndexes const& branchListIndexes() const;

    EventToProcessBlockIndexes const& eventToProcessBlockIndexes() const;

    Provenance const& getProvenance(ProductID const& pid) const;
    StableProvenance const& getStableProvenance(ProductID const& pid) const;

    BasicHandle getByProductID(ProductID const& oid) const;

    void put(BranchDescription const& bd,
             std::unique_ptr<WrapperBase> edp,
             ProductProvenance const& productProvenance) const;

    void put(ProductResolverIndex index, std::unique_ptr<WrapperBase> edp, ParentageID productProvenance) const;

    void putOnRead(BranchDescription const& bd,
                   std::unique_ptr<WrapperBase> edp,
                   std::optional<ProductProvenance> productProvenance) const;

    WrapperBase const* getIt(ProductID const& pid) const override;
    std::optional<std::tuple<WrapperBase const*, unsigned int>> getThinnedProduct(ProductID const& pid,
                                                                                  unsigned int key) const override;
    void getThinnedProducts(ProductID const& pid,
                            std::vector<WrapperBase const*>& foundContainers,
                            std::vector<unsigned int>& keys) const override;
    OptionalThinnedKey getThinnedKeyFrom(ProductID const& parent,
                                         unsigned int key,
                                         ProductID const& thinned) const override;

    ProductID branchIDToProductID(BranchID const& bid) const;

    void mergeProvenanceRetrievers(EventPrincipal& other) {
      provRetrieverPtr_->mergeProvenanceRetrievers(other.provRetrieverPtr());
    }

    using Base::getProvenance;
    using Base::getStableProvenance;

    unsigned int processBlockIndex(std::string const& processName) const override;

  private:
    EventPrincipal(std::shared_ptr<ProductRegistry const> reg,
                   std::vector<std::shared_ptr<ProductResolverBase>>&& resolvers,
                   std::shared_ptr<BranchIDListHelper const> branchIDListHelper,
                   std::shared_ptr<ThinnedAssociationsHelper const> thinnedAssociationsHelper,
                   ProcessConfiguration const& pc,
                   HistoryAppender* historyAppender,
                   unsigned int streamIndex,
                   ProcessBlockHelperBase const* processBlockHelper);

    BranchID pidToBid(ProductID const& pid) const;

    edm::ThinnedAssociation const* getThinnedAssociation(edm::BranchID const& branchID) const;

    unsigned int transitionIndex_() const override;
    void changedIndexes_() final;

    std::shared_ptr<ProductProvenanceRetriever const> provRetrieverPtr() const {
      return get_underlying_safe(provRetrieverPtr_);
    }
    std::shared_ptr<ProductProvenanceRetriever>& provRetrieverPtr() { return get_underlying_safe(provRetrieverPtr_); }

    bool wasBranchListIndexesChangedFromInput(BranchListIndexes const&) const;
    void updateBranchListIndexes(BranchListIndexes&&);
    void commonFillEventPrincipal(EventAuxiliary const& aux,
                                  ProcessHistory const* processHistory,
                                  DelayedReader* reader);

  private:
    EventAuxiliary aux_;

    edm::propagate_const<LuminosityBlockPrincipal*> luminosityBlockPrincipal_;

    // Pointer to the 'retriever' that will get provenance information from the persistent store.
    edm::propagate_const<std::shared_ptr<ProductProvenanceRetriever>> provRetrieverPtr_;

    EventSelectionIDVector eventSelectionIDs_;

    std::shared_ptr<BranchIDListHelper const> branchIDListHelper_;
    ProcessBlockHelperBase const* processBlockHelper_;
    std::shared_ptr<ThinnedAssociationsHelper const> thinnedAssociationsHelper_;

    BranchListIndexes branchListIndexes_;

    EventToProcessBlockIndexes eventToProcessBlockIndexes_;

    std::vector<ProcessIndex> branchListIndexToProcessIndex_;

    StreamID streamID_;
  };

  inline bool isSameEvent(EventPrincipal const& a, EventPrincipal const& b) { return isSameEvent(a.aux(), b.aux()); }
}  // namespace edm
#endif
