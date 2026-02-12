//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// count_min_sketch.cpp
//
// Identification: src/primer/count_min_sketch.cpp
//
// Copyright (c) 2015-2025, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "primer/count_min_sketch.h"

#include <queue>
#include <stdexcept>
#include <string>
#include <utility>

namespace bustub {

/**
 * Constructor for the count-min sketch.
 *
 * @param width The width of the sketch matrix.
 * @param depth The depth of the sketch matrix.
 * @throws std::invalid_argument if width or depth are zero.
 */
template <typename KeyType>
CountMinSketch<KeyType>::CountMinSketch(uint32_t width, uint32_t depth) : width_(width), depth_(depth) {
  if (width_ == 0 || depth_ == 0) {
    throw std::invalid_argument("Width or depth cannot be 0.");
  }
  this->buckets.resize(depth_);
  for (size_t i = 0; i < depth_; i++) {
    buckets[i] = std::vector<std::atomic<uint32_t>>(width);
    for (size_t j = 0; j < width_; j++) {
      buckets[i][j].store(0, std::memory_order_relaxed);
    }
  }
  /** @spring2026 PLEASE DO NOT MODIFY THE FOLLOWING */
  // Initialize seeded hash functions
  hash_functions_.reserve(depth_);
  for (size_t i = 0; i < depth_; i++) {
    hash_functions_.push_back(this->HashFunction(i));
  }
}

template <typename KeyType>
CountMinSketch<KeyType>::CountMinSketch(CountMinSketch &&other) noexcept
                                        : width_(other.width_)
                                        , depth_(other.depth_)
                                        , buckets(std::move(other.buckets))
                                        , hash_functions_(std::move(other.hash_functions_)) {
}

template <typename KeyType>
auto CountMinSketch<KeyType>::operator=(CountMinSketch &&other) noexcept -> CountMinSketch & {
  if (this != &other) {
    width_ = other.width_;
    depth_ = other.depth_;
    buckets = std::move(other.buckets);
    hash_functions_ = std::move(other.hash_functions_);
  }
  return *this;
}

template <typename KeyType>
void CountMinSketch<KeyType>::Insert(const KeyType &item) {
  // mutex_.lock();
  for (size_t i = 0; i < depth_; i++) {
    buckets[i][hash_functions_[i](item) % width_].fetch_add(1, std::memory_order_relaxed);
  }
  // mutex_.unlock();
}

template <typename KeyType>
void CountMinSketch<KeyType>::Merge(const CountMinSketch<KeyType> &other) {
  if (width_ != other.width_ || depth_ != other.depth_) {
    throw std::invalid_argument("Incompatible CountMinSketch dimensions for merge.");
  }
  for (size_t i = 0; i < depth_; i++) {
    for (size_t j = 0; j < width_; j++) {
      buckets[i][j].fetch_add(other.buckets[i][j].load(std::memory_order_relaxed),
        std::memory_order_relaxed);
    }
  }
  /** @TODO(student) Implement this function! */
}

template <typename KeyType>
auto CountMinSketch<KeyType>::Count(const KeyType &item) const -> uint32_t {
  uint32_t min = buckets[0][hash_functions_[0](item) % width_];
  for (size_t i = 1; i < depth_; i++) {
    uint32_t idx = hash_functions_[i](item) % width_;
    uint32_t val = buckets[i][idx].load(std::memory_order_relaxed);
    if (val < min) {
      min = val;
    }
  }
  return min;
}

template <typename KeyType>
void CountMinSketch<KeyType>::Clear() {
  /** @TODO(student) Implement this function! */
  for (size_t i = 0; i < depth_; i++) {
    for (size_t j = 0; j < width_; j++) {
      buckets[i][j].store(0, std::memory_order_relaxed);
    }
  }
}

template <typename KeyType>
auto CountMinSketch<KeyType>::TopK(uint16_t k, const std::vector<KeyType> &candidates)
    -> std::vector<std::pair<KeyType, uint32_t>> {
  /** @TODO(student) Implement this function! */
  using Pair = std::pair<KeyType, uint32_t>;
  if (k == 0 || candidates.empty()) {
    return {};
  }
  auto cmp = [](const Pair &a, const Pair &b) {
    return a.second > b.second;
  };
  std::priority_queue<Pair, std::vector<Pair>, decltype(cmp)> min_heap(cmp);
  for (auto &candidate: candidates) {
    uint32_t count = Count(candidate);
    if (min_heap.size() < k) {
      min_heap.emplace(candidate, count);
    } else if (min_heap.top().second < count) {
      min_heap.pop();
      min_heap.emplace(candidate, count);
    }
  }
  std::vector<Pair> results;
  results.reserve(min_heap.size());
  while (!min_heap.empty()) {
    results.push_back(min_heap.top());
    min_heap.pop();
  }
  std::sort(results.begin(), results.end(),
    [](const Pair &a, const Pair &b) {
      return a.second > b.second;
    });
  return results;
}

// Explicit instantiations for all types used in tests
template class CountMinSketch<std::string>;
template class CountMinSketch<int64_t>;  // For int64_t tests
template class CountMinSketch<int>;      // This covers both int and int32_t
}  // namespace bustub
