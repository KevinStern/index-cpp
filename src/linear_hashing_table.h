/* Copyright (c) 2014 Kevin L. Stern
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

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

// Non-portable Linux-specific headers for file truncation.
#include <unistd.h>
#include <sys/types.h>

#include "bitarray.h"
#include "buffer.h"

// Linear hashing is a hash table technique which is designed to incur the cost of expansion and
// contraction in small increments, as items are added and removed, instead of with less frequent
// but more costly operations that copy the entire underlying data store.
// The incremental expansion/contraction operations of LinearHashingTable make it ideal for database
// applications. In particular, since the load factor of the hash table tends to hover near
// the thresholds during growth and contraction phases, the table tends not to reserve significant
// amounts of empty space within the underlying storage. In addition, the table never requires the
// copy of the entire underlying data store, an operation that would severely limit the size of the
// table.
//
// Data in LinearHashingTable are organized into buckets, each of which holds a configurable number
// of records (KeyType/ValueType pairs). Primary buckets are indexed by hash code and overflow
// buckets form a chain off of a primary bucket as buckets within the hash category are filled. When
// the load factor of the hash table exceeds a configurable threshold, a new primary bucket is added
// to the table and the records within a bucket chain are split among itself and the new bucket
// chain. When the load factor of the hash table falls below a configurable threshold, the most
// recently added bucket is removed and its records are added to the bucket chain that was split
// when it was newly added.
//
// At construction time, LinearHashingTable accepts an instance of BucketManager.
// TransientBucketManager gives in-memory storage behavior while FileBucketManager gives filesystem
// storage behavior.
template <typename KeyType, typename ValueType, typename HashEvaluator = std::hash<KeyType>,
          typename EqualityEvaluator = std::equal_to<KeyType>>
class LinearHashingTable {
private:
  typedef std::pair<KeyType, ValueType> RecordType;

  // A struct holding information about the structure of the table.
  //
  // The struct is maintained by the bucket managers.
  struct TableMetadata {
    // The capacity of each bucket. Buckets can hold a number of records not exceeding the capacity.
    size_t bucket_capacity;
    // The minimum number of buckets composing the table.
    size_t min_bucket_count;
    // n is the current base bucket count (a power of two) and p points to the next bucket to be
    // split. Note that when all buckets have been split, n is doubled and p is reset to zero.
    size_t n, p;
    // The minimum and maximum load factors.
    double min_load_factor, max_load_factor;
    // The number of records across all buckets in the table.
    size_t size;
    // The number of primary buckets.
    size_t primary_bucket_count;

    // Constructor which initializes with default minimum bucket count and load factors, and the
    // specified bucket capacity.
    TableMetadata(size_t bucket_capacity)
        : bucket_capacity(bucket_capacity), min_bucket_count(32), n(32), p(0),
          min_load_factor(0.75), max_load_factor(0.95), size(0), primary_bucket_count(0) {}

    // Constructor which initializes with default load factors, and the specified bucket capacity
    // and minimum bucket count. Note that min_bucket_count must be a power of two.
    TableMetadata(size_t bucket_capacity, size_t min_bucket_count)
        : bucket_capacity(bucket_capacity), min_bucket_count(min_bucket_count), n(min_bucket_count),
          p(0), min_load_factor(0.75), max_load_factor(0.95), size(0), primary_bucket_count(0) {
      if ((min_bucket_count & (min_bucket_count - 1)) != 0) {
        throw std::invalid_argument("min_bucket_count not a power of two");
      }
    }

    // Constructor which initializes with the specified load factors, minimum bucket count and
    // bucket capacity. Note that min_bucket_count must be a power of two.
    TableMetadata(size_t bucket_capacity, size_t min_bucket_count, double min_load_factor,
                  double max_load_factor)
        : bucket_capacity(bucket_capacity), min_bucket_count(min_bucket_count), n(min_bucket_count),
          p(0), min_load_factor(min_load_factor), max_load_factor(max_load_factor), size(0),
          primary_bucket_count(0) {
      if ((min_bucket_count & (min_bucket_count - 1)) != 0) {
        throw std::invalid_argument("min_bucket_count not a power of two");
      }
    }
  };

public:
  // BucketManager creates, flushes and loads buckets, and maintains TableMetadata for
  // LinearHashingTable.
  class BucketManager {
  public:
    // Bucket is the atom of data storage behind LinearHashingTable, it holds records as well as a
    // pointer to its overflow bucket.
    class Bucket {
    public:
      // Convenience class which locks flushing for a bucket until it falls out of scope.
      class FlushLock {
      public:
        FlushLock(Bucket* bucket) : bucket_(bucket) { bucket->lock_flush(); }
        ~FlushLock() { bucket_->unlock_flush(); }
      private:
        Bucket* bucket_;
      };

      // index: A user-specified index value carried along by this bucket. FileBucketManager, e.g.,
      //        uses this to store the index at which this bucket is stored in the underlying
      //        storage, while BucketManager (its in-memory counterpart) ignores this value.
      Bucket(BucketManager* bucket_manager, size_t index, bool is_overflow, bool is_dirty)
          : bucket_manager_(bucket_manager), index_(index), is_overflow_(is_overflow),
            is_dirty_(is_dirty), flush_lock_count_(0) {}

      ~Bucket() {}

      // Chain a new overflow bucket off of this bucket. If this bucket already has an overflow
      // bucket, an error is generated.
      // Will result in a bucket load, if this bucket is unloaded.
      void create_overflow() {
        FlushLock lock(this);
        bucket_manager_->load(this, records_, overflow_);
        if (overflow_ != nullptr) {
          throw std::logic_error("overflow bucket already exists");
        }
        overflow_.reset(bucket_manager_->new_overflow_bucket());
        set_dirty(true);
      }

      // Access the data held by this bucket.
      // Will result in a bucket load, if this bucket is unloaded.
      const std::vector<RecordType>& data() {
        bucket_manager_->load(this, records_, overflow_);
        return records_;
      }

      // Access the data element at index i.
      // Will result in a bucket load, if this bucket is unloaded.
      RecordType& data(size_t i) {
        bucket_manager_->load(this, records_, overflow_);
        return records_[i];
      }

      // Get the identifying index value provided at construction time.
      size_t index() const {
        return index_;
      }

      // Will result in a bucket load, if this bucket is unloaded.
      void insert(const KeyType& key, const ValueType& value) {
        bucket_manager_->load(this, records_, overflow_);
        records_.push_back(RecordType(key, value));
        set_dirty(true);
      }

      // Whether or not this bucket has un-flushed changes.
      bool is_dirty() const {
        return is_dirty_;
      }

      // Whether or not this bucket is an overflow bucket.
      bool is_overflow() const {
        return is_overflow_;
      }

      // Whether or not this bucket may be flushed from the cache. If the number of calls to
      // lock_flush() exceeds the number of calls to unlock_flush, then flush is locked;
      // otherwise, it is not.
      bool is_flush_locked() const {
        return flush_lock_count_ > 0;
      }

      // Will result in a bucket load, if this bucket is unloaded.
      bool is_full() {
        return size() == bucket_manager_->metadata_.bucket_capacity;
      }

      // Control whether or not this bucket may be flushed from the cache. If the number of calls
      // to lock_flush() exceeds the number of calls to unlock_flush, then flush is locked;
      // otherwise, it is not.
      // FlushLock is a convenient mechanism for locking a bucket as it will unlock the Bucket when
      // it falls out of scope.
      void lock_flush() {
        ++flush_lock_count_;
      }

      // Get the overflow bucket chained off of this Bucket; nullptr if none.
      Bucket* overflow() const {
        return overflow_.get();
      }

      // Remove the element at index i from this bucket.
      // Will result in a bucket load, if this bucket is unloaded.
      void remove(size_t i) {
        bucket_manager_->load(this, records_, overflow_);
        records_.erase(records_.begin() + i);
        set_dirty(true);
      }

      // Set whether or not this bucket has unflushed changes.
      void set_dirty(bool is_dirty) {
        is_dirty_ = is_dirty;
      }

      // Get the number of elements in this bucket.
      // Will result in a bucket load, if this bucket is unloaded.
      size_t size() {
        bucket_manager_->load(this, records_, overflow_);
        return records_.size();
      }

      // Unload this bucket. Clears all records stored in memory.
      void unload() {
        std::vector<RecordType> empty;
        records_.swap(empty);
      }

      // Control whether or not this bucket may be flushed from the cache. If the number of calls
      // to lock_flush() exceeds the number of calls to unlock_flush, then flush is locked;
      // otherwise, it is not.
      // FlushLock is a convenient mechanism for locking a bucket as it will unlock the Bucket when
      // it falls out of scope.
      void unlock_flush() {
        if (flush_lock_count_ == 0) {
          throw std::logic_error("flush not locked");
        }
        --flush_lock_count_;
      }

    private:
      BucketManager* bucket_manager_;
      size_t index_;
      bool is_overflow_;
      bool is_dirty_;
      size_t flush_lock_count_;
      std::vector<RecordType> records_;
      std::unique_ptr<Bucket> overflow_;
    };

    // bucket_capacity: The maximum number of records that each bucket can hold.
    BucketManager(size_t bucket_capacity) : metadata_(bucket_capacity) {}

    // bucket_capacity: The maximum number of records that each bucket can hold.
    // min_bucket_count: The minimum number of buckets composing the table.
    BucketManager(size_t bucket_capacity, size_t min_bucket_count)
        : metadata_(bucket_capacity, min_bucket_count) {}

    // bucket_capacity: The maximum number of records that each bucket can hold.
    // min_bucket_count: The minimum number of buckets composing the table.
    // min_load_factor: The minimum load factor below which capacity will be removed.
    // max_load_factor: The maximum load factor above which capacity will be added.
    BucketManager(size_t bucket_capacity, size_t min_bucket_count, double min_load_factor,
                  double max_load_factor)
        : metadata_(bucket_capacity, min_bucket_count, min_load_factor, max_load_factor) {}

    virtual ~BucketManager() {}

    // Remove bucket from any record-keeping storage. Call with a bucket when removing it from the
    // table.
    virtual void discard(Bucket* bucket) {}

    // Flush state to storage.
    virtual void flush() = 0;

    // Initialize table metadata.
    virtual void initialize() = 0;

    // Load bucket by populating its records and its overflow bucket, if necessary.
    virtual void load(Bucket* bucket, std::vector<RecordType>& records,
                      std::unique_ptr<Bucket>& overflow) = 0;

    // Create a new, empty bucket at the next available index.
    virtual Bucket* new_overflow_bucket() = 0;

    // Create a new, empty bucket at the specified index.
    virtual Bucket* new_primary_bucket(size_t index) = 0;

    // Create an unloaded bucket for index.
    virtual Bucket* unloaded_bucket(size_t index, bool is_overflow) = 0;

  protected:
    TableMetadata metadata_;

    friend class LinearHashingTable;
  };
  typedef typename BucketManager::Bucket Bucket;
  typedef typename Bucket::FlushLock FlushLock;

  // A BucketManager that gives an in-memory only table. It does not persist buckets.
  class TransientBucketManager : public BucketManager {
  public:
    // bucket_capacity: The maximum number of records that each bucket can hold.
    TransientBucketManager(size_t bucket_capacity) : BucketManager(bucket_capacity) {}

    // bucket_capacity: The maximum number of records that each bucket can hold.
    // min_bucket_count: The minimum number of buckets composing the table.
    TransientBucketManager(size_t bucket_capacity, size_t min_bucket_count)
        : BucketManager(bucket_capacity, min_bucket_count) {}

    // bucket_capacity: The maximum number of records that each bucket can hold.
    // min_bucket_count: The minimum number of buckets composing the table.
    // min_load_factor: The minimum load factor below which capacity will be removed.
    // max_load_factor: The maximum load factor above which capacity will be added.
    TransientBucketManager(size_t bucket_capacity, size_t min_bucket_count, double min_load_factor,
                  double max_load_factor)
        : BucketManager(bucket_capacity, min_bucket_count, min_load_factor, max_load_factor) {}

    virtual void discard(Bucket* bucket) override {}

    virtual void flush() override {}

    virtual void initialize() override {}

    virtual void load(Bucket* bucket, std::vector<RecordType>& records,
                      std::unique_ptr<Bucket>& overflow) override {}

    virtual Bucket* new_overflow_bucket() override {
      return new Bucket(this, 0, true, false);
    }

    virtual Bucket* new_primary_bucket(size_t index) override {
      return new Bucket(this, index, false, false);
    }

    virtual Bucket* unloaded_bucket(size_t index, bool is_overflow) override {
      return new Bucket(this, index, is_overflow, false);
    }
  };

  // A BucketManager that stores buckets across two files, a primary bucket file and an overflow
  // bucket file. The primary bucket file holds both the primary buckets and the table metadata
  // while the overflow bucket file holds only the overflow buckets. Neither file need exist before
  // the manager is constructed.
  class FileBucketManager : public BucketManager {
  public:
    // bucket_capacity: The maximum number of records that each bucket can hold.
    // primary_filename: The name of the file holding the primary buckets.
    // overflow_filename: The name of the file holding the overflow buckets.
    // serialize: A routine for serializing records.
    // deserialize: A routine for deserializing records.
    // encoded_record_size: The size of a serialized record, in bytes.
    // endian: The endianness with which to store integral values.
    FileBucketManager(size_t bucket_capacity, const char* primary_filename,
        const char* overflow_filename,
        std::function<void(const KeyType, const ValueType, Buffer&)> serialize,
        std::function<void(Buffer&, KeyType&, ValueType&)> deserialize, size_t encoded_record_size,
        Buffer::Endian endian, size_t bucket_cache_size)
        : BucketManager(bucket_capacity), serialize_(serialize), deserialize_(deserialize),
          encoded_record_size_(encoded_record_size), endian_(endian),
          bucket_cache_size_(bucket_cache_size) {
      primary_file_ = open_file(primary_filename);
      overflow_file_ = open_file(overflow_filename);
    }

    // bucket_capacity: The maximum number of records that each bucket can hold.
    // min_bucket_count: The minimum number of buckets composing the table.
    // primary_filename: The name of the file holding the primary buckets.
    // overflow_filename: The name of the file holding the overflow buckets.
    // serialize: A routine for serializing records.
    // deserialize: A routine for deserializing records.
    // encoded_record_size: The size of a serialized record, in bytes.
    // endian: The endianness with which to store integral values.
    FileBucketManager(size_t bucket_capacity, size_t min_bucket_count, const char* primary_filename,
        const char* overflow_filename,
        std::function<void(const KeyType, const ValueType, Buffer&)> serialize,
        std::function<void(Buffer&, KeyType&, ValueType&)> deserialize, size_t encoded_record_size,
        Buffer::Endian endian, size_t bucket_cache_size)
        : BucketManager(bucket_capacity, min_bucket_count), serialize_(serialize),
          deserialize_(deserialize), encoded_record_size_(encoded_record_size), endian_(endian),
          bucket_cache_size_(bucket_cache_size) {
      primary_file_ = open_file(primary_filename);
      overflow_file_ = open_file(overflow_filename);
    }

    // bucket_capacity: The maximum number of records that each bucket can hold.
    // min_bucket_count: The minimum number of buckets composing the table.
    // min_load_factor: The minimum load factor below which capacity will be removed.
    // max_load_factor: The maximum load factor above which capacity will be added.
    // primary_filename: The name of the file holding the primary buckets.
    // overflow_filename: The name of the file holding the overflow buckets.
    // serialize: A routine for serializing records.
    // deserialize: A routine for deserializing records.
    // encoded_record_size: The size of a serialized record, in bytes.
    // endian: The endianness with which to store integral values.
    FileBucketManager(size_t bucket_capacity, size_t min_bucket_count, double min_load_factor,
        double max_load_factor, const char* primary_filename, const char* overflow_filename,
        std::function<void(const KeyType, const ValueType, Buffer&)> serialize,
        std::function<void(Buffer&, KeyType&, ValueType&)> deserialize, size_t encoded_record_size,
        Buffer::Endian endian, size_t bucket_cache_size)
        : BucketManager(bucket_capacity, min_bucket_count, min_load_factor, max_load_factor),
          serialize_(serialize), deserialize_(deserialize),
          encoded_record_size_(encoded_record_size), endian_(endian),
          bucket_cache_size_(bucket_cache_size) {
      primary_file_ = open_file(primary_filename);
      overflow_file_ = open_file(overflow_filename);
    }

    virtual ~FileBucketManager() {
      fclose(primary_file_);
      fclose(overflow_file_);
    }

    virtual void discard(Bucket* bucket) override {
      loaded_.erase(bucket);
      if (bucket->is_overflow()) {
        overflow_bucket_status_.unset(bucket->index());
      }
    }

    // Flush all dirty buckets to file and unload all buckets from the cache that are not locked for
    // flush.
    virtual void flush() override {
      bool changes_flushed = false;
      for (auto i = loaded_.begin(); i != loaded_.end();) {
        Bucket* bucket = *i;
        if (bucket->is_flush_locked()) {
          ++i;
          continue;
        }
        if (bucket->is_dirty()) {
          FILE* file = bucket->is_overflow() ? overflow_file_ : primary_file_;
          bucket->set_dirty(false);
          changes_flushed = true;
          const std::vector<RecordType>& records = bucket->data();
          Buffer buffer(endian_, BUCKET_METADATA_SIZE + records.size() * encoded_record_size_);
          Bucket* overflow = bucket->overflow();
          if (overflow != nullptr) {
            buffer.put_uint64(overflow->index());
          } else {
            buffer.put_uint64(std::numeric_limits<uint64_t>::max());
          }
          buffer.put_uint64(records.size());
          for (const auto& next : records) {
            serialize_(next.first, next.second, buffer);
          }
          buffer.switch_mode();
          if (fseek(file, to_file_offset(bucket->index()), SEEK_SET) != 0) {
            throw std::runtime_error("seek error");
          }
          if (fwrite(buffer.data(), 1, buffer.remaining(), file) != buffer.remaining()) {
            throw std::runtime_error("write error");
          }
        }
        bucket->unload();
        i = loaded_.erase(i);
      }
      if (changes_flushed) {
        flush_table_metadata();
      }
    }

    // Reads the table metadata from the primary storage file.
    virtual void initialize() override {
      if (fseek(primary_file_, 0, SEEK_END) != 0) {
        throw std::runtime_error("seek error");
      }
      size_t primary_size = ftell(primary_file_);
      if (primary_size == 0) {
        return;
      }
      if (fseek(primary_file_, -TABLE_METADATA_SIZE, SEEK_END) != 0) {
        throw std::runtime_error("seek error");
      }

      Buffer metabuf(endian_, TABLE_METADATA_SIZE);
      metabuf.inc_index(fread(metabuf.data(), 1, metabuf.remaining(), primary_file_));
      metabuf.switch_mode();
      size_t overflow_bucket_size = metabuf.get_uint64();
      size_t overflow_bucket_count = Bitarray::byte_count(overflow_bucket_size);
      BucketManager::metadata_.bucket_capacity = metabuf.get_uint64();
      BucketManager::metadata_.min_bucket_count = metabuf.get_uint64();
      BucketManager::metadata_.n = metabuf.get_uint64();
      BucketManager::metadata_.p = metabuf.get_uint64();
      BucketManager::metadata_.min_load_factor = metabuf.get_double();
      BucketManager::metadata_.max_load_factor = metabuf.get_double();
      BucketManager::metadata_.size = metabuf.get_uint64();
      BucketManager::metadata_.primary_bucket_count = metabuf.get_uint64();
      if (fseek(primary_file_, -(TABLE_METADATA_SIZE + overflow_bucket_count), SEEK_END) != 0) {
        throw std::runtime_error("seek error");
      }
      Buffer buf(endian_, overflow_bucket_count);
      buf.inc_index(fread(buf.data(), 1, buf.remaining(), primary_file_));
      buf.switch_mode();
      std::vector<uint8_t> overflow_bucket_status_data;
      while (!buf.exhausted()) {
        overflow_bucket_status_data.push_back(buf.get_ubyte());
      }
      overflow_bucket_status_.set_bytes(overflow_bucket_status_data, overflow_bucket_size);
    }

    virtual void load(Bucket* bucket, std::vector<RecordType>& records,
                      std::unique_ptr<Bucket>& overflow) override {
      if (!loaded_.insert(bucket).second) {
        return;
      }
      FILE* file = bucket->is_overflow() ? overflow_file_ : primary_file_;
      Buffer metabuf(endian_, BUCKET_METADATA_SIZE);

      if (fseek(file, to_file_offset(bucket->index()), SEEK_SET) != 0) {
        throw std::runtime_error("seek error");
      }
      metabuf.inc_index(fread(metabuf.data(), 1, metabuf.remaining(), file));
      metabuf.switch_mode();
      size_t overflow_index = metabuf.get_uint64();
      size_t count = metabuf.get_uint64();

      Buffer buffer(endian_, count * encoded_record_size_);
      buffer.inc_index(fread(buffer.data(), 1, buffer.remaining(), file));
      buffer.switch_mode();
      records.resize(count);
      for (int i = 0; i < count; ++i) {
        deserialize_(buffer, records[i].first, records[i].second);
      }
      if (bucket->overflow() == nullptr) {
        overflow.reset(
            overflow_index == std::numeric_limits<uint64_t>::max() ?
                nullptr :
                unloaded_bucket(overflow_index, true));
      }
      if (loaded_.size() > bucket_cache_size_) {
        FlushLock lock(bucket);
        flush();
      }
    }

    virtual Bucket* new_overflow_bucket() override {
      size_t i = 0;
      for (; i < overflow_bucket_status_.size(); ++i) {
        if (!overflow_bucket_status_.is_set(i)) {
          break;
        }
      }
      if (i == overflow_bucket_status_.size()) {
        overflow_bucket_status_.resize(i + 1);
      }
      overflow_bucket_status_.set(i);

      Bucket* bucket = new Bucket(this, i, true, true);
      loaded_.insert(bucket);
      if (loaded_.size() > bucket_cache_size_) {
        FlushLock lock(bucket);
        flush();
      }
      return bucket;
    }

    virtual Bucket* new_primary_bucket(size_t index) override {
      Bucket* bucket = new Bucket(this, index, false, true);
      loaded_.insert(bucket);
      if (loaded_.size() > bucket_cache_size_) {
        FlushLock lock(bucket);
        flush();
      }
      return bucket;
    }

    virtual Bucket* unloaded_bucket(size_t index, bool is_overflow) override {
      return new Bucket(this, index, is_overflow, false);
    }

  private:
    static const int TABLE_METADATA_SIZE = 72;
    static const int BUCKET_METADATA_SIZE = 16;
    static const int FREE = 0;
    static const int OCCUPIED = 1;

    static FILE* open_file(const char* filename) {
      FILE* file = fopen(filename, "rb+");
      if (file == nullptr) {
        file = fopen(filename, "wb+");
        if (file == nullptr) {
          throw std::runtime_error("file open error");
        }
      }
      if (setvbuf(file, nullptr, _IONBF, 0) != 0) {
        throw std::runtime_error("io buffer error");
      }
      return file;
    }

    // Flush the table metadata to the primary storage file.
    void flush_table_metadata() {
      // Crop unused overflow buckets.
      size_t overflow_bucket_count;
      for (overflow_bucket_count = overflow_bucket_status_.size();
           overflow_bucket_count > 0 && !overflow_bucket_status_.is_set(overflow_bucket_count - 1);
           --overflow_bucket_count);
      overflow_bucket_status_.resize(overflow_bucket_count);

      const std::vector<uint8_t>& overflow_bucket_status_data = overflow_bucket_status_.bytes();
      Buffer buffer(endian_, TABLE_METADATA_SIZE + overflow_bucket_status_data.size());

      buffer.put_bytes(&overflow_bucket_status_data[0], overflow_bucket_status_data.size());

      buffer.put_uint64(overflow_bucket_status_.size());
      buffer.put_uint64(BucketManager::metadata_.bucket_capacity);
      buffer.put_uint64(BucketManager::metadata_.min_bucket_count);
      buffer.put_uint64(BucketManager::metadata_.n);
      buffer.put_uint64(BucketManager::metadata_.p);
      buffer.put_double(BucketManager::metadata_.min_load_factor);
      buffer.put_double(BucketManager::metadata_.max_load_factor);
      buffer.put_uint64(BucketManager::metadata_.size);
      buffer.put_uint64(BucketManager::metadata_.primary_bucket_count);
      buffer.switch_mode();
      if (fseek(primary_file_,
                to_file_offset(BucketManager::metadata_.primary_bucket_count), SEEK_SET) != 0) {
        throw std::runtime_error("seek error");
      }
      if (fwrite(buffer.data(), 1, buffer.remaining(), primary_file_) != buffer.remaining()) {
        throw std::runtime_error("write error");
      }
      // These are non-portable, but there is no portable way to in-place truncate a file.
      if (ftruncate(fileno(primary_file_), ftell(primary_file_)) != 0) {
        throw std::runtime_error("truncate error");
      }
      if (ftruncate(fileno(overflow_file_), to_file_offset(overflow_bucket_status_.size())) != 0) {
        throw std::runtime_error("truncate error");
      }
    }

    // Returns the offset within the storage file of the bucket at index bucket_index.
    size_t to_file_offset(size_t bucket_index) const {
      size_t result = bucket_index *
          (BUCKET_METADATA_SIZE + BucketManager::metadata_.bucket_capacity * encoded_record_size_);
      return result;
    }

    std::function<void(const KeyType, const ValueType, Buffer&)> serialize_;
    std::function<void(Buffer&, KeyType&, ValueType&)> deserialize_;
    size_t encoded_record_size_;
    Buffer::Endian endian_;
    size_t bucket_cache_size_;
    Bitarray overflow_bucket_status_;
    std::unordered_set<Bucket*> loaded_;

    FILE* primary_file_;
    FILE* overflow_file_;
  };

  // bucket_manager: The BucketManager to use for constructing buckets of records; client retains
  //                 ownership of bucket_manager.
  LinearHashingTable(BucketManager* bucket_manager)
      : bucket_manager_(bucket_manager), metadata_(bucket_manager->metadata_) {
    initialize();
  }

  // bucket_manager: The BucketManager to use for constructing buckets of records; client retains
  //                 ownership of bucket_manager.
  LinearHashingTable(BucketManager* bucket_manager, HashEvaluator hash, EqualityEvaluator equals)
      : bucket_manager_(bucket_manager), hash_(hash), equals_(equals),
        metadata_(bucket_manager->metadata_) {
    initialize();
  }

  ~LinearHashingTable() {
    bucket_manager_->flush();
  }

  // Insert a mapping from key to value into the table, if no mapping from key already exists.
  // Returns true if a mapping from key to value was inserted into the table, false otherwise.
  bool insert(const KeyType& key, const ValueType& value) {
    size_t bucket_index = calculate_bucket_chain(key);
    Bucket* dummy_bucket;
    size_t dummy_index;
    if (find_internal(key, bucket_index, dummy_bucket, dummy_index)) {
      return false;
    }
    insert_internal(key, value, bucket_index);
    ++metadata_.size;
    rebalance_split();
    return true;
  }

  // Remove the mapping from key from the table and populate value with its target, if such a
  // mapping exists.
  // Returns true if a mapping from key was removed from the table, false otherwise.
  bool remove(const KeyType& key, ValueType& value) {
    Bucket* bucket;
    size_t index;
    if (find_internal(key, calculate_bucket_chain(key), bucket, index)) {
      value = bucket->data(index).second;
      bucket->remove(index);
      --metadata_.size;
      rebalance_merge();
      return true;
    }
    return false;
  }

  // Test whether or not the table contains a mapping from key.
  bool contains(const KeyType& key) {
    Bucket* dummy_bucket;
    size_t dummy_index;
    return find_internal(key, calculate_bucket_chain(key), dummy_bucket, dummy_index);
  }

  // Find the mapping from key and populate value with its target.
  // Returns true if a mapping from key was found, false otherwise.
  bool find(const KeyType& key, ValueType& value) {
    Bucket* bucket;
    size_t index;
    if (find_internal(key, calculate_bucket_chain(key), bucket, index)) {
      value = bucket->data(index).second;
      return true;
    }
    return false;
  }

  // Get the number of records in the table.
  size_t size() {
    return metadata_.size;
  }

private:
  BucketManager* bucket_manager_;
  HashEvaluator hash_;
  EqualityEvaluator equals_;
  TableMetadata& metadata_;
  std::vector<std::unique_ptr<Bucket>> bucket_chains_;

  // Calculate the bucket chain to which key belongs.
  size_t calculate_bucket_chain(const KeyType& key) const {
    size_t digest = hash_(key);
    size_t bucket = digest & (metadata_.n - 1);
    if (bucket < metadata_.p) {
      bucket = digest & ((metadata_.n << 1) - 1);
    }
    return bucket;
  }

  // Compute the current table load. The table load is defined as the ratio of the current size to
  // the sum of the capacities of the primary buckets across all bucket chains.
  double compute_load() const {
    return static_cast<double>(metadata_.size) /
        (bucket_manager_->metadata_.bucket_capacity * bucket_chains_.size());
  }

  // Search for key in the primary bucket chain at bucket_index. If found, populate source with the
  // containing bucket, index with key's position within source, and return true; otherwise, return
  // false.
  bool find_internal(const KeyType& key, size_t bucket_index, Bucket*& source, size_t& index) {
    source = bucket_chains_[bucket_index].get();
    while (source != nullptr) {
      const std::vector<RecordType>& records = source->data();
      auto i = std::find_if(records.begin(), records.end(),
                            [&](const RecordType& r)->bool{return equals_(r.first, key);});
      if (i != records.end()) {
        index = i - records.begin();
        return true;
      }
      source = source->overflow();
    }
    return false;
  }

  // Expand the number of primary buckets by adding a bucket to the end of the sequence of bucket
  // chains.
  virtual void grow() {
    ++metadata_.primary_bucket_count;
    bucket_chains_.push_back(std::unique_ptr<Bucket>(
        bucket_manager_->new_primary_bucket(bucket_chains_.size())));
  }

  // Initialize the bucket manager and construct the bucket chains.
  void initialize() {
    bucket_manager_->initialize();
    if (metadata_.primary_bucket_count > 0) {
      for (int i = 0; i < metadata_.primary_bucket_count; ++i) {
        // Create the primary bucket.
        bucket_chains_.push_back(
            std::unique_ptr<Bucket>(bucket_manager_->unloaded_bucket(i, false)));
      }
    } else {
      for (int i = 0; i < metadata_.n; ++i) grow();
    }
  }

  // Insert key+value into the primary bucket chain at bucket_index.
  void insert_internal(const KeyType& key, const ValueType& value, size_t bucket_index) {
    Bucket* source = bucket_chains_[bucket_index].get();
    while (source->is_full()) {
      if (source->overflow() == nullptr) {
        source->create_overflow();
      }
      source = source->overflow();
    }
    source->insert(key, value);
  }

  // Checks whether or not the current load exceeds the maximum load factor; if so, resolves by
  // splitting the bucket chain at index p_, repeating until the current load falls below the
  // maximum load factor.
  void rebalance_split() {
    while (compute_load() > metadata_.max_load_factor) {
      size_t split_index = metadata_.p;
      // Load the bucket to be split and lock flush.
      FlushLock l1(bucket_chains_[split_index].get());
      bucket_chains_[split_index]->data();
      std::unique_ptr<Bucket> target(split(split_index));
      // Lock flush to the replacement of the split bucket.
      FlushLock l2(bucket_chains_[split_index].get());
      if (++metadata_.p == metadata_.n) {
        metadata_.p = 0;
        metadata_.n <<= 1;
      }
      Bucket* bucket = target.get();
      while (bucket != nullptr) {
        FlushLock lock(bucket);
        const std::vector<RecordType>& records = bucket->data();
        for (const auto& next : records) {
          insert_internal(next.first, next.second, calculate_bucket_chain(next.first));
        }
        bucket_manager_->discard(bucket);
        bucket = bucket->overflow();
      }
    }
  }

  // Checks whether or not the current load is under the minimum load factor; if so, resolves by
  // merging the bucket chains resulting from the most recent split, repeating until the current
  // load rises above the minimum load factor unless the number of bucket chains would fall below
  // the minimum size of the table.
  void rebalance_merge() {
    while (bucket_chains_.size() > metadata_.min_bucket_count &&
           compute_load() < metadata_.min_load_factor) {
      std::unique_ptr<Bucket> target(shrink());
      if (metadata_.p == 0) {
        metadata_.n >>= 1;
        metadata_.p = metadata_.n - 1;
      } else {
        --metadata_.p;
      }
      Bucket* bucket = target.get();
      while (bucket != nullptr) {
        FlushLock lock(bucket);
        const std::vector<RecordType>& records = bucket->data();
        for (const auto& next : records) {
          insert_internal(next.first, next.second, calculate_bucket_chain(next.first));
        }
        bucket_manager_->discard(bucket);
        bucket = bucket->overflow();
      }
    }
  }

  // Contract the number of primary buckets by removing a bucket from the end of the sequence of
  // bucket chains.
  // Returns the bucket that was removed to the caller, who then owns the Bucket.
  virtual Bucket* shrink() {
    --metadata_.primary_bucket_count;
    Bucket* result = bucket_chains_.back().release();
    bucket_chains_.pop_back();
    return result;
  }

  // Swap the bucket to be split with a replacement bucket and add its companion.
  // Returns the bucket that was split to the caller, who then owns the Bucket.
  Bucket* split(size_t bucket_index) {
    Bucket* bucket = bucket_manager_->new_primary_bucket(bucket_chains_[bucket_index]->index());
    FlushLock lock(bucket);
    Bucket* result = bucket_chains_[bucket_index].release();
    bucket_chains_[bucket_index].reset(bucket);
    grow();
    return result;
  }
};
