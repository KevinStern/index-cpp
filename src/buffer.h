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

#include <boost/detail/endian.hpp>

#include <algorithm>
#include <cstring>
#include <istream>
#include <stdexcept>

// A fixed size buffer of bytes.
// 
// Exposes setters and getters for encoding and decoding integral and floating point numbers in
// addition to those for single and multi byte collections, strings and streams.
//
// General maneuvering through the buffer involves manipulating three elements:
// 1) The <def>index</def> is the current read/write position within the buffer; a get or put
//    operation will occur at the current index so long as there is room between the current index
//    and the stop.
// 2) The <def>stop</def> restricts the index of the buffer for get and put type operations, it can
//    be manipulated on-demand, but may never exceed the capacity of the buffer.
// 3) The <def>capacity</def> is the total amount of memory allocated for the buffer.
class Buffer {
public:
  enum Endian {
    LITTLE, BIG
  };

  // capacity: The capacity of the buffer.
  Buffer(uint64_t capacity) : index_(0), stop_(capacity), capacity_(capacity),
      bytes_(new uint8_t[capacity]), endian_mismatch_(false) {}

  // order: The endianness with which to write integral types.
  // capacity: The capacity of the buffer.
  Buffer(Endian order, uint64_t capacity) :
      index_(0), stop_(capacity), capacity_(capacity), bytes_(new uint8_t[capacity]),
#ifdef BOOST_BIG_ENDIAN
      endian_mismatch_(order == LITTLE)
#else
      endian_mismatch_(order == BIG)
#endif
  {}

  ~Buffer() {
    delete [] bytes_;
  }

  // Fill the buffer with data from the stream.
  void fill_from(std::istream* stream) {
    stream->read(reinterpret_cast<char*>(bytes_ + index_), remaining());
    index_ += stream->gcount();
  }

  // Copy count elements of chars to the buffer.
  void put_chars(const char* chars, uint64_t count) {
    put_bytes(reinterpret_cast<const uint8_t*>(chars), count);
  }

  // Copy count elements of bytes to the buffer.
  void put_bytes(const int8_t* bytes, uint64_t count) {
    put_bytes(reinterpret_cast<const uint8_t*>(bytes), count);
  }

  // Copy count elements of bytes to the buffer.
  void put_bytes(const uint8_t* bytes, uint64_t count) {
    if (index_ > stop_ - count) {
      throw std::overflow_error("overflow");
    }
    memcpy(bytes_ + index_, bytes, count);
    index_ += count;
  }

  // Copy val to the buffer.
  void put_byte(int8_t val) {
    internal_put_integer8(val);
  }

  // Copy val to the buffer.
  void put_ubyte(uint8_t val) {
    internal_put_integer8(val);
  }

  // Copy val to the buffer.
  void put_int32(int32_t val) {
    internal_put_integer32(val);
  }

  // Copy val to the buffer.
  void put_uint32(uint32_t val) {
    internal_put_integer32(val);
  }

  // Copy val to the buffer.
  void put_int64(int64_t val) {
    internal_put_integer64(val);
  }

  // Copy val to the buffer.
  void put_uint64(uint64_t val) {
    internal_put_integer64(val);
  }

  // Copy val to the buffer.
  void put_float(float val) {
    if (index_ > stop_ - 4) {
      throw std::overflow_error("overflow");
    }
    uint8_t* p = reinterpret_cast<uint8_t*>(&val);
    memcpy(bytes_ + index_, p, 4);
    index_ += 4;
  }

  // Copy val to the buffer.
  void put_double(double val) {
    if (index_ > stop_ - 8) {
      throw std::overflow_error("overflow");
    }
    uint8_t* p = reinterpret_cast<uint8_t*>(&val);
    memcpy(bytes_ + index_, p, 8);
    index_ += 8;
  }

  // Drain the buffer to stream.
  void drain_to(std::ostream& stream) {
    uint64_t count = remaining();
    stream.write(reinterpret_cast<char*>(bytes_ + index_), count);
    index_ += count;
  }

  // Drain up to count bytes from the buffer to stream; if count > remaining() then only remaining()
  // bytes will be drained.
  uint64_t drain_to(std::ostream& stream, uint64_t count) {
    uint64_t rem = remaining();
    if (rem < count) {
      count = rem;
    }
    stream.write(reinterpret_cast<char*>(bytes_ + index_), count);
    index_ += count;
    return count;
  }

  // Copy count elements of the buffer to target.
  void get_chars(char* target, uint64_t count) {
    get_bytes(reinterpret_cast<uint8_t*>(target), count);
  }

  // Copy count elements of the buffer to target.
  void get_bytes(int8_t* target, uint64_t count) {
    get_bytes(reinterpret_cast<uint8_t*>(target), count);
  }

  // Copy count elements of the buffer to target.
  void get_bytes(uint8_t* target, uint64_t count) {
    if (index_ > stop_ - count) {
      throw std::underflow_error("underflow");
    }
    memcpy(target, bytes_ + index_, count);
    index_ += count;
  }

  // Get the next byte from the buffer.
  int8_t get_byte() {
    return internal_get_integer8<int8_t>();
  }

  // Get the next unsigned byte from the buffer.
  uint8_t get_ubyte() {
    return internal_get_integer8<uint8_t>();
  }

  // Get the next int32 from the buffer.
  int32_t get_int32() {
    return internal_get_integer32<int32_t>();
  }

  // Get the next uint32 from the buffer.
  uint32_t get_uint32() {
    return internal_get_integer32<uint32_t>();
  }

  // Get the next int64 from the buffer.
  int64_t get_int64() {
    return internal_get_integer64<int64_t>();
  }

  // Get the next uint64 from the buffer.
  uint64_t get_uint64() {
    return internal_get_integer64<uint64_t>();
  }

  // Get the next float from the buffer.
  float get_float() {
    if (index_ > stop_ - 4) {
      throw std::underflow_error("underflow");
    }
    float result = *reinterpret_cast<float*>(bytes_ + index_);
    index_ += 4;
    return result;
  }

  // Get the next double from the buffer.
  double get_double() {
    if (index_ > stop_ - 8) {
      throw std::underflow_error("underflow");
    }
    double result = *reinterpret_cast<double*>(bytes_ + index_);
    index_ += 8;
    return result;
  }

  // Get the index.
  uint64_t index() {
    return index_;
  }

  // Set the index.
  //
  // index: The new index; must not exceed the stop.
  void set_index(uint64_t index) {
    if (index > stop_) {
      throw std::invalid_argument("index exceeds stop");
    }
    index_ = index;
  }

  // Increment the index by delta.
  void inc_index(uint64_t delta) {
    index_ += delta;
  }

  // Get the stop.
  uint64_t stop() {
    return stop_;
  }

  // Set the stop.
  //
  // stop: The new stop; must not exceed the capacity of the buffer.
  void set_stop(uint64_t stop) {
    if (stop > capacity_) {
      throw std::invalid_argument("stop exceeds capacity");
    }
    stop_ = stop;
  }

  // Get the capacity of the buffer.
  uint64_t capacity() {
    return capacity_;
  }

  // Returns true if the index is equal to the stop; false otherwise.
  bool exhausted() {
    return index_ == stop_;
  }

  // Set the index to the stop so that remaining() returns false.
  void exhaust() {
    index_ = stop_;
  }

  // Get the number of bytes between the index and the stop.
  uint64_t remaining() {
    return stop_ - index_;
  }

  // Set the stop to the index and the index to zero.
  //
  // Useful, e.g., to switch from write mode to read mode.
  void switch_mode() {
    set_stop(index_);
    set_index(0);
  }

  // Set the stop to the capacity of the buffer and the index to zero.
  void reset() {
    set_stop(capacity_);
    set_index(0);
  }

  // Get a direct reference to the internal data store at the current index.
  uint8_t* data() {
    return bytes_ + index_;
  }

private:
  uint64_t index_, stop_, capacity_;
  uint8_t* bytes_;
  bool endian_mismatch_;

  // Disable default
  Buffer(const Buffer& buffer) {}

  // Disable default
  Buffer& operator=(const Buffer& buffer) {return *this;}

  template <typename T>
  inline T internal_get_integer8() {
    if (index_ > stop_ - 1) {
      throw std::underflow_error("underflow");
    }
    T result = *reinterpret_cast<T*>(bytes_ + index_);
    index_ += 1;
    return result;
  }

  template <typename T>
  inline T internal_get_integer32() {
    if (index_ > stop_ - 4) {
      throw std::underflow_error("underflow");
    }
    T result = *reinterpret_cast<T*>(bytes_ + index_);
    if (endian_mismatch_) {
      swap32(reinterpret_cast<uint8_t*>(&result));
    }
    index_ += 4;
    return result;
  }

  template <typename T>
  inline T internal_get_integer64() {
    if (index_ > stop_ - 8) {
      throw std::underflow_error("underflow");
    }
    T result = *reinterpret_cast<T*>(bytes_ + index_);
    if (endian_mismatch_) {
      swap64(reinterpret_cast<uint8_t*>(&result));
    }
    index_ += 8;
    return result;
  }

  template <typename T>
  inline void internal_put_integer8(T byte) {
    if (index_ > stop_ - 1) {
      throw std::overflow_error("overflow");
    }
    uint8_t* p = reinterpret_cast<uint8_t*>(&byte);
    bytes_[index_] = *p;
    index_ += 1;
  }

  template <typename T>
  inline void internal_put_integer32(T quad) {
    if (index_ > stop_ - 4) {
      throw std::overflow_error("overflow");
    }
    uint8_t* p = reinterpret_cast<uint8_t*>(&quad);
    if (endian_mismatch_) {
      swap32(p);
    }
    memcpy(bytes_ + index_, p, 4);
    index_ += 4;
  }

  template <typename T>
  inline void internal_put_integer64(T oct) {
    if (index_ > stop_ - 8) {
      throw std::overflow_error("overflow");
    }
    uint8_t* p = reinterpret_cast<uint8_t*>(&oct);
    if (endian_mismatch_) {
      swap64(p);
    }
    memcpy(bytes_ + index_, p, 8);
    index_ += 8;
  }

  inline void swap(uint8_t* p1, uint8_t* p2) {
    uint8_t temp = *p1;
    *p1 = *p2;
    *p2 = temp;
  }

  inline void swap32(uint8_t* p) {
    swap(p, p + 3);
    swap(p + 1, p + 2);
  }

  inline void swap64(uint8_t* p) {
    swap(p, p + 7);
    swap(p + 1, p + 6);
    swap(p + 2, p + 5);
    swap(p + 3, p + 4);
  }
};
