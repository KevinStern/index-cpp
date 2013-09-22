/* Copyright (c) 2013 Kevin L. Stern
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#pragma once

#include <chrono>
#include <cstdint>
#include <functional>

#include "buffer.h"
#include "geoprimitives.h"
#include "record.h"

// GeoIndex is an index for GeoPoints within CoordinateSystem that provides fast range queries.
template <typename RecordType, typename CoordinateSystem>
class GeoIndex {
private:
  static const int DEFAULT_PREFIX_REFINEMENT_STOP_SIZE = 8;

public:
  class DataSource {
  public:
    // Given function lt that returns true for all records at indices less than i and false for all
    // records at indices greater than or equal to i, jump to the record at index i.
    virtual void jump_to(std::function<bool(const RecordType&)> lt) const = 0;

    // Initialize the next record from the data source; this record is accessed via record().
    // Returns false if the data source is exhausted; true otherwise.
    virtual bool next() const = 0;

    // Get the record read at the latest call to next().
    virtual const RecordType& record() const = 0;
  };

  // DataSource that pulls from a file that contains a sequence of serialized RecordTypes, ordered
  // by their GeoHash value.
  class FileDataSource : public DataSource {
  public:
    // filename: The geo index file.
    // deserialize: Function for parsing a RecordType from a Buffer.
    // encoded_record_size: The number of bytes that a serialized instance of RecordType consumes.
    // endian: The endianness of integral types encoded in the record file.
    // buffer_size_records: The size of the read buffer in number of records.
    FileDataSource(const char* filename, std::function<void(Buffer&, RecordType&)> deserialize,
        uint32_t encoded_record_size, Buffer::Endian endian, uint32_t buffer_size_records)
        : reader_(filename, deserialize, encoded_record_size, endian, buffer_size_records) {}

    virtual void jump_to(std::function<bool(const RecordType&)> lt) const override {
      reader_.jump_to(lt);
    }

    virtual bool next() const override {
      return reader_.read(record_);
    }

    virtual const RecordType& record() const override {
      return record_;
    }

  private:
    mutable RecordReader<RecordType> reader_;
    mutable RecordType record_;
  };

  // DataSource that pulls from a vector that contains a sequence of RecordTypes, ordered by their
  // GeoHash value.
  class VectorDataSource : public DataSource {
  public:
    // records: A vector of records, odered by their GeoHash value.
    VectorDataSource(const std::vector<RecordType>& records) : records_(records), index_(-1) {}

    virtual void jump_to(std::function<bool(const RecordType&)> lt) const override {
      index_ = lower_bound(records_.begin(), records_.end(), records_[0],
          [&](const RecordType& r1, const RecordType& r2)->bool{return lt(r1);}) - records_.begin()
          - 1;
    }

    virtual bool next() const override {
      return ++index_ < records_.size();
    }

    virtual const RecordType& record() const override {
      return records_[index_];
    }

  private:
    mutable std::vector<RecordType> records_;
    mutable int64_t index_;
  };

  struct GeoRecord {
    GeoHash<CoordinateSystem> hash;
    GeoPoint<CoordinateSystem> point;

    GeoRecord(const GeoPoint<CoordinateSystem>& point)
        : hash(Hasher().hash(point)), point(point) {}
  };

  // Statistics to give visibility into the geo-index query process.
  struct QueryStats {
    uint64_t record_read_count;
    uint64_t prefix_search_count;
    uint64_t time_us;

    QueryStats() : record_read_count(0), prefix_search_count(0), time_us(0) {}
  };

  // Order instances of RecordType by their GeoHash value.
  // begin/end: Iterator giving the RecordTypes to order.
  // geopoint: Function for producing GeoPoints from instances of RecordType.
  // target: Storage for holding the ordered RecordTypes; may be the same container that produced
  //         the iterator.
  template <typename Iterator>
  static void order(Iterator begin, Iterator end,
                    std::function<GeoPoint<CoordinateSystem>(const RecordType&)> geopoint,
                    std::vector<RecordType>& target) {
    std::vector<CombinedRecord> c;
    while (begin != end) {
      const RecordType& record = *begin;
      GeoRecord geo(geopoint(record));
      c.push_back(CombinedRecord(geo, record));
      ++begin;
    }
    target.clear();
    std::sort(c.begin(), c.end(),
        [](const CombinedRecord& r1, const CombinedRecord& r2)->bool{
          return r1.first.hash.value() < r2.first.hash.value();
        });
    target.reserve(c.size());
    for (auto& next : c) {
      target.push_back(next.second);
    }
  }

  // source: The DataSource from which to pull index records; client retains ownership over source.
  // geopoint: Function for producing GeoPoints from instances of RecordType.
  // prefix_refinement_stop_size: The prefix count at which to stop refining query prefixes.
  GeoIndex(const DataSource* source,
           std::function<GeoPoint<CoordinateSystem>(const RecordType&)> geopoint,
           int prefix_refinement_stop_size = DEFAULT_PREFIX_REFINEMENT_STOP_SIZE)
      : source_(source), geopoint_(geopoint),
        prefix_refinement_stop_size_(prefix_refinement_stop_size) {}

  // Query the geo index for records that fall within bounds.
  QueryStats query(const GeoRectangle<CoordinateSystem>& bounds,
      std::function<void(const RecordType&)> callback) const {
    time_t time = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();

    QueryStats stats;
    std::vector<GeoPrefix<CoordinateSystem>> prefixes;
    GeoPrefix<CoordinateSystem>::search_prefixes(bounds, prefixes, prefix_refinement_stop_size_);
    stats.prefix_search_count = prefixes.size();
    for (const auto& prefix : prefixes) {
      source_->jump_to([&](const RecordType& record)->bool{
        GeoRecord geo(geopoint_(record));
        return geo.hash.value() < prefix.value();
      });
      RecordType record;
      while (source_->next()) {
        const RecordType& record = source_->record();
        GeoRecord geo(geopoint_(record));
        if (!prefix.applies_to(geo.hash.value())) {
          break;
        }
        ++stats.record_read_count;
        if (bounds.contains(geo.point)) {
          callback(record);
        }
      }
    }

    stats.time_us = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count() - time;

    return stats;
  }

private:
  typedef std::pair<GeoRecord, RecordType> CombinedRecord;

  static const GeoHasher& Hasher() {
    static GeoHasher hasher;
    return hasher;
  }

  const DataSource* source_;
  std::function<GeoPoint<CoordinateSystem>(const RecordType&)> geopoint_;
  int prefix_refinement_stop_size_;
};
