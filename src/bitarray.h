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

#include <cstdint>
#include <istream>
#include <stdexcept>
#include <vector>

// A dynamic array of bits.
class Bitarray {
public:
  // Calculate the number of bytes required to store size bits.
  static size_t byte_count(size_t size) {
    return (size >> 3) + 1;
  }

  Bitarray() : size_(0) {}

  // Get a read only view of the underlying data.
  // Bytes are organized so that index 0 is at bit 0 of byte 0, while index n is at bit
  // n & ((1 << 3) - 1) of byte n >> 3.
  const std::vector<uint8_t>& bytes() const {
    return bytes_;
  }

  // Test whether or not the bit at index is set.
  bool is_set(size_t index) const {
    if (index >= size_) {
      throw std::invalid_argument("index out of range");
    }
    size_t byte = index >> 3;
    uint8_t bit = index & ((1 << 3) - 1);
    return (bytes_[byte] & (1 << bit));
  }

  // Resize the array to size bits. If the array must expand to accommodate the size, new unset bits
  // are added.
  void resize(size_t size) {
    bytes_.resize(byte_count(size), 0);
    size_ = size;
  }

  // Set the bit at index.
  void set(size_t index) {
    if (index >= size_) {
      throw std::invalid_argument("index out of range");
    }
    size_t byte = index >> 3;
    uint8_t bit = index & ((1 << 3) - 1);
    bytes_[byte] |= (1 << bit);
  }

  // Set the size and the underlying data.
  // Bytes are organized so that index 0 is at bit 0 of byte 0, while index n is at bit
  // n & ((1 << 3) - 1) of byte n >> 3.
  void set_bytes(std::vector<uint8_t> bytes, size_t size) {
    bytes_.swap(bytes);
    size_ = size;
  }

  // Get the number of bits in the array.
  size_t size() const {
    return size_;
  }

  // Unset the bit at index.
  void unset(size_t index) {
    if (index >= size_) {
      throw std::invalid_argument("index out of range");
    }
    size_t byte = index >> 3;
    uint8_t bit = index & ((1 << 3) - 1);
    bytes_[byte] &= ~(1 << bit);
  }

private:
  size_t size_;
  std::vector<uint8_t> bytes_;

  // For pretty printing.
  friend std::ostream& operator<<(std::ostream& out, const Bitarray& bitarray) {
    out << "[";
    if (bitarray.size_ > 0) {
      for (int i = 0; ; ++i) {
        out << bitarray.is_set(i);
        if (i < bitarray.size_ - 1) {
          out << ",";
        } else {
          break;
        }
      }
    }
    out << "]";
    return out;
  }
};
