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

#include <cstdint>
#include <cstdio>
#include <stdexcept>

#include "buffer.h"

// A class for reading a file containing an ordered list of records.
//
// A client specifies a filename and gives a record parser along with other information about the
// encoding of instances of RecordType and uses this class for record access. Two access mechanisms
// are provided: indexed access and binary search. Indexed access consumes constant time per
// retrieval and binary search consumes logarithmic time in the number of records per retrieval.
template <typename RecordType>
class RecordReader {
public:
  // Construct a reader to draw from the record file indicated by file.
  // 
  // file: The path to the record file.
  // parser: A function for parsing a record from a Buffer.
  // encoded_record_size: The size in bytes of an encoded record.
  // endian: The endianness of integral types encoded in the record file.
  // buffer_size_records: The size of the read buffer in number of records.
  RecordReader(const char* file, std::function<void(Buffer&, RecordType&)> parser,
               uint32_t encoded_record_size, Buffer::Endian endian,
               uint32_t buffer_size_records) :
               parser_(parser), encoded_record_size_(encoded_record_size),
               buffer_(endian, encoded_record_size * buffer_size_records), index_(0) {
    file_ = fopen(file, "rb");
    if (file_ == nullptr) {
      throw std::runtime_error("file open error");
    }
    if (fseek(file_, 0, SEEK_END) != 0) {
      throw std::runtime_error("seek error");
    }
    size_t file_size;
    if ((file_size = ftell(file_)) < 0) {
      throw std::runtime_error("position error");
    }
    record_count_ = file_size / encoded_record_size_;
    if (fseek(file_, 0, SEEK_SET) != 0) {
      throw std::runtime_error("seek error");
    }
    if (setvbuf(file_, nullptr, _IONBF, 0) != 0) {
      throw std::runtime_error("io buffer error");
    }
    buffer_.exhaust();
  }

  ~RecordReader() {
    if (fclose(file_) != 0) {
      throw std::runtime_error("file close error");
    }
  }

  // Get the index of the next record to be returned by read.
  uint64_t index() const {
    return index_;
  }

  // Jump to the record at the specified index.
  void jump(uint64_t index) {
    if (index > record_count_) {
      throw std::invalid_argument("index exceeds record count");
    }
    internal_jump(index);
    buffer_.exhaust();
  }

  // Given function lt that returns true for all records at indices less than i and false for all
  // records at indices greater than or equal to i, jump to index i.
  // 
  // This function operates in time logarithmic in the number of records.
  void jump_to(std::function<bool(const RecordType&)> lt) {
    // Configure the buffer to read a single record.
    buffer_.set_stop(encoded_record_size_);

    RecordType record;
    int64_t lo = 0;
    int64_t hi = record_count_;
    while (lo < hi) {
      int64_t mid = (lo + hi) >> 1;
      buffer_.set_index(0);
      internal_jump(mid);
      internal_fill_buffer();
      buffer_.set_index(0);
      parser_(buffer_, record);
      if (lt(record)) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    internal_jump(hi);
    buffer_.exhaust();
  }

  // Parse the next record into record.
  bool read(RecordType& record) {
    if (buffer_.remaining() == 0) {
      buffer_.reset();
      internal_fill_buffer();
      buffer_.switch_mode();
    }
    if (buffer_.remaining() == 0) {
      return false;
    }
    parser_(buffer_, record);
    ++index_;
    return true;
  }

  // Get the record count.
  uint64_t record_count() {
    return record_count_;
  }

private:
  FILE* file_;
  std::function<void(Buffer&, RecordType&)> parser_;
  const uint32_t encoded_record_size_;

  Buffer buffer_;
  uint64_t record_count_;

  uint64_t index_;

  void internal_fill_buffer() {
    size_t count = fread(buffer_.data(), 1, buffer_.remaining(), file_);
    buffer_.inc_index(count);
  }

  void internal_jump(uint64_t index) {
    index_ = index;
    if (fseek(file_, index * encoded_record_size_, SEEK_SET) != 0) {
      throw std::runtime_error("seek error");
    }
  }
};
