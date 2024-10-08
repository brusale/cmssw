#ifndef FWCore_Framework_global_producerAbilityToImplementor_h
#define FWCore_Framework_global_producerAbilityToImplementor_h
// -*- C++ -*-
//
// Package:     FWCore/Framework
// File  :     producerAbilityToImplementor
//
/**\file  producerAbilityToImplementor.h "FWCore/Framework/interface/global/producerAbilityToImplementor.h"

 Description: Class used to pair a module Ability to the actual base class used to implement that ability

 Usage:
    <usage>

*/
//
// Original Author:  Chris Jones
//         Created:  Thu, 18 Jul 2013 11:51:33 GMT
//

// system include files

// user include files
#include "FWCore/Framework/interface/moduleAbilities.h"
#include "FWCore/Framework/interface/global/implementors.h"
#include "FWCore/Framework/interface/global/EDProducerBase.h"

// forward declarations
namespace edm {
  namespace global {
    namespace producer {
      template <typename T>
      struct AbilityToImplementor;

      template <typename C>
      struct AbilityToImplementor<edm::StreamCache<C>> {
        using Type = edm::global::impl::StreamCacheHolder<edm::global::EDProducerBase, C>;
      };

      template <typename... Cs>
      struct AbilityToImplementor<edm::InputProcessBlockCache<Cs...>> {
        using Type = edm::global::impl::InputProcessBlockCacheHolder<edm::global::EDProducerBase, Cs...>;
      };

      template <typename C>
      struct AbilityToImplementor<edm::RunCache<C>> {
        using Type = edm::global::impl::RunCacheHolder<edm::global::EDProducerBase, C>;
      };

      template <typename C>
      struct AbilityToImplementor<edm::RunSummaryCache<C>> {
        using Type = edm::global::impl::RunSummaryCacheHolder<edm::global::EDProducerBase, C>;
      };

      template <typename C>
      struct AbilityToImplementor<edm::LuminosityBlockCache<C>> {
        using Type = edm::global::impl::LuminosityBlockCacheHolder<edm::global::EDProducerBase, C>;
      };

      template <typename C>
      struct AbilityToImplementor<edm::LuminosityBlockSummaryCache<C>> {
        using Type = edm::global::impl::LuminosityBlockSummaryCacheHolder<edm::global::EDProducerBase, C>;
      };

      template <>
      struct AbilityToImplementor<edm::WatchProcessBlock> {
        using Type = edm::global::impl::WatchProcessBlock<edm::global::EDProducerBase>;
      };

      template <>
      struct AbilityToImplementor<edm::BeginProcessBlockProducer> {
        using Type = edm::global::impl::BeginProcessBlockProducer<edm::global::EDProducerBase>;
      };

      template <>
      struct AbilityToImplementor<edm::EndProcessBlockProducer> {
        using Type = edm::global::impl::EndProcessBlockProducer<edm::global::EDProducerBase>;
      };

      template <>
      struct AbilityToImplementor<edm::BeginRunProducer> {
        using Type = edm::global::impl::BeginRunProducer<edm::global::EDProducerBase>;
      };

      template <>
      struct AbilityToImplementor<edm::EndRunProducer> {
        using Type = edm::global::impl::EndRunProducer<edm::global::EDProducerBase>;
      };

      template <>
      struct AbilityToImplementor<edm::BeginLuminosityBlockProducer> {
        using Type = edm::global::impl::BeginLuminosityBlockProducer<edm::global::EDProducerBase>;
      };

      template <>
      struct AbilityToImplementor<edm::EndLuminosityBlockProducer> {
        using Type = edm::global::impl::EndLuminosityBlockProducer<edm::global::EDProducerBase>;
      };

      template <>
      struct AbilityToImplementor<edm::ExternalWork> {
        using Type = edm::global::impl::ExternalWork<edm::global::EDProducerBase>;
      };

      template <>
      struct AbilityToImplementor<edm::Transformer> {
        using Type = edm::global::impl::Transformer<edm::global::EDProducerBase>;
      };

      template <>
      struct AbilityToImplementor<edm::Accumulator> {
        using Type = edm::global::impl::Accumulator<edm::global::EDProducerBase>;
      };

      template <bool, bool, typename T>
      struct SpecializeAbilityToImplementor {
        using Type = typename AbilityToImplementor<T>::Type;
      };

      template <bool B, typename C>
      struct SpecializeAbilityToImplementor<true, B, edm::RunSummaryCache<C>> {
        using Type = typename edm::global::impl::EndRunSummaryProducer<edm::global::EDProducerBase, C>;
      };

      template <bool B>
      struct SpecializeAbilityToImplementor<true, B, edm::EndRunProducer> {
        using Type = edm::global::impl::EmptyType;
      };

      template <bool B, typename C>
      struct SpecializeAbilityToImplementor<B, true, edm::LuminosityBlockSummaryCache<C>> {
        using Type = typename edm::global::impl::EndLuminosityBlockSummaryProducer<edm::global::EDProducerBase, C>;
      };

      template <bool B>
      struct SpecializeAbilityToImplementor<B, true, edm::EndLuminosityBlockProducer> {
        using Type = edm::global::impl::EmptyType;
      };
    }  // namespace producer
  }  // namespace global
}  // namespace edm

#endif
